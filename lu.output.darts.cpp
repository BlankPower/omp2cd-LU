#include "lu.output.darts.h"
using namespace darts;
using namespace std;
std::mutex TP9992mutex;
static boolean  flag [13] ;
static void verify(double xcr[5], double xce[5], double xci, char *class_is, boolean *verified);
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
int main(int argc, char **argv) {
getOMPNumThreads();
getOMPSchedulePolicy();
getTPLoopThresholds();
getNumTPs();
affin  = new ThreadAffinity(ompNumThreads/NUMTPS - 1, NUMTPS, COMPACT, getDARTSTPPolicy(), getDARTSMCPolicy());affinMaskRes = affin->generateMask ();
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
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP192>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (int*)&((nthreads))));
}
ssor();
printf("ssor complete...\n");
error();
printf("error complete...\n");
pintgr();
printf("pintgr complete...\n");
verify(rsdnm, errnm, frc, &class_is, &verified);
mflops = (double)itmax * (1984.77 * (double)nx0 * (double)ny0 * (double)nz0 - 10923.299999999999 * (((double)(nx0 + ny0 + nz0) / 3.) * ((double)(nx0 + ny0 + nz0) / 3.)) + 27770.900000000001 * (double)(nx0 + ny0 + nz0) / 3. - 144010.) / (maxtime * 1.0E+6);
printf("verify complete...\n");
c_print_results("LU", class_is, nx0, ny0, nz0, itmax, nthreads, maxtime, mflops, "          floating point", verified, "3.0 structured", "09 Sep 2022", "gcc", "gcc", "-lm -fopenmp", "-I../common", "-O3 -fopenmp", "(none)", "(none)");
}
/*Function: domain, ID: 3*/
static void domain() {
/*domain:3*/
/*CompoundStmt:2183*/
nx = nx0;
ny = ny0;
nz = nz0;
if (nx < 4 || ny < 4 || nz < 4) {
    printf("     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 4 THEY ARE CURRENTLY%3d%3d%3d\n", nx, ny, nz);
    exit(1);
}
if (nx > 12 || ny > 12 || nz > 12) {
    printf("     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%4d%4d%4d\n", nx, ny, nz);
    exit(1);
}
ist = 1;
iend = nx - 2;
jst = 1;
jend = ny - 2;
}
/*Function: erhs, ID: 4*/
static void erhs() {
/*erhs:4*/
/*CompoundStmt:2203*/
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP2204>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
}
}
/*Function: error, ID: 5*/
static void error() {
/*error:5*/
/*CompoundStmt:4730*/
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
static void exact(int i, int j, int k, double u000ijk[5]) {
/*exact:6*/
/*CompoundStmt:4794*/
int m;
double xi, eta, zeta;
xi = ((double)i) / (nx0 - 1);
eta = ((double)j) / (ny0 - 1);
zeta = ((double)k) / (nz - 1);
for (m = 0; m < 5; m++) {
    u000ijk[m] = ce[m][0] + ce[m][1] * xi + ce[m][2] * eta + ce[m][3] * zeta + ce[m][4] * xi * xi + ce[m][5] * eta * eta + ce[m][6] * zeta * zeta + ce[m][7] * xi * xi * xi + ce[m][8] * eta * eta * eta + ce[m][9] * zeta * zeta * zeta + ce[m][10] * xi * xi * xi * xi + ce[m][11] * eta * eta * eta * eta + ce[m][12] * zeta * zeta * zeta * zeta;
}
}
/*Function: l2norm, ID: 9*/
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double sum[5]) {
/*l2norm:9*/
/*CompoundStmt:9878*/
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP9879>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (int*)&((iend)), (int*)&((ist)), (int*)&((jend)), (int*)&((jst)), (int*)&((nx0)), (int*)&((ny0)), (int*)&((nz0)), (double **)&((sum))));
}
}
/*Function: pintgr, ID: 10*/
static void pintgr() {
/*pintgr:10*/
/*CompoundStmt:10021*/
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
        phi1[i][j] = 0.40000000000000002 * (u[i][j][k][4] - 0.5 * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2])) + ((u[i][j][k][3]) * (u[i][j][k][3]))) / u[i][j][k][0]);
        k = ki2;
        phi2[i][j] = 0.40000000000000002 * (u[i][j][k][4] - 0.5 * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2])) + ((u[i][j][k][3]) * (u[i][j][k][3]))) / u[i][j][k][0]);
    }
}
frc1 = 0.;
for (i = ibeg; i <= ifin1; i++) {
    for (j = jbeg; j <= jfin1; j++) {
        frc1 = frc1 + (phi1[i][j] + phi1[i + 1][j] + phi1[i][j + 1] + phi1[i + 1][j + 1] + phi2[i][j] + phi2[i + 1][j] + phi2[i][j + 1] + phi2[i + 1][j + 1]);
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
            phi1[i][k] = 0.40000000000000002 * (u[i][jbeg][k][4] - 0.5 * (((u[i][jbeg][k][1]) * (u[i][jbeg][k][1])) + ((u[i][jbeg][k][2]) * (u[i][jbeg][k][2])) + ((u[i][jbeg][k][3]) * (u[i][jbeg][k][3]))) / u[i][jbeg][k][0]);
        }
    }
}
jglob = jfin;
if (jglob == ji2) {
    for (i = ibeg; i <= ifin; i++) {
        iglob = i;
        for (k = ki1; k <= ki2; k++) {
            phi2[i][k] = 0.40000000000000002 * (u[i][jfin][k][4] - 0.5 * (((u[i][jfin][k][1]) * (u[i][jfin][k][1])) + ((u[i][jfin][k][2]) * (u[i][jfin][k][2])) + ((u[i][jfin][k][3]) * (u[i][jfin][k][3]))) / u[i][jfin][k][0]);
        }
    }
}
frc2 = 0.;
for (i = ibeg; i <= ifin1; i++) {
    for (k = ki1; k <= ki2 - 1; k++) {
        frc2 = frc2 + (phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1] + phi2[i][k] + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1]);
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
            phi1[j][k] = 0.40000000000000002 * (u[ibeg][j][k][4] - 0.5 * (((u[ibeg][j][k][1]) * (u[ibeg][j][k][1])) + ((u[ibeg][j][k][2]) * (u[ibeg][j][k][2])) + ((u[ibeg][j][k][3]) * (u[ibeg][j][k][3]))) / u[ibeg][j][k][0]);
        }
    }
}
iglob = ifin;
if (iglob == ii2) {
    for (j = jbeg; j <= jfin; j++) {
        jglob = j;
        for (k = ki1; k <= ki2; k++) {
            phi2[j][k] = 0.40000000000000002 * (u[ifin][j][k][4] - 0.5 * (((u[ifin][j][k][1]) * (u[ifin][j][k][1])) + ((u[ifin][j][k][2]) * (u[ifin][j][k][2])) + ((u[ifin][j][k][3]) * (u[ifin][j][k][3]))) / u[ifin][j][k][0]);
        }
    }
}
frc3 = 0.;
for (j = jbeg; j <= jfin1; j++) {
    for (k = ki1; k <= ki2 - 1; k++) {
        frc3 = frc3 + (phi1[j][k] + phi1[j + 1][k] + phi1[j][k + 1] + phi1[j + 1][k + 1] + phi2[j][k] + phi2[j + 1][k] + phi2[j][k + 1] + phi2[j + 1][k + 1]);
    }
}
frc3 = deta * dzeta * frc3;
frc = 0.25 * (frc1 + frc2 + frc3);
}
/*Function: read_input, ID: 11*/
static void read_input() {
/*read_input:11*/
/*CompoundStmt:10651*/
FILE *fp;
printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - LU Benchmark\n\n");
fp = fopen("inputlu.data", "r");
if (fp != ((void *)0)) {
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
    printf("     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n");
    exit(1);
}
if (nx0 > 12 || ny0 > 12 || nz0 > 12) {
    printf("     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
    exit(1);
}
printf(" Size: %3dx%3dx%3d\n", nx0, ny0, nz0);
printf(" Iterations: %3d\n", itmax);
}
/*Function: rhs, ID: 12*/
static void rhs() {
/*rhs:12*/
/*CompoundStmt:10795*/
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP10796>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
}
}
/*Function: setbv, ID: 13*/
static void setbv() {
/*setbv:13*/
/*CompoundStmt:13196*/
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP13197>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
}
}
/*Function: setcoeff, ID: 14*/
static void setcoeff() {
/*setcoeff:14*/
/*CompoundStmt:13424*/
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
static void setiv() {
/*setiv:15*/
/*CompoundStmt:13763*/
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP13764>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
}
}
/*Function: ssor, ID: 16*/
static void ssor() {
/*ssor:16*/
/*CompoundStmt:13885*/
int i, j, k, m;
int istep;
double tmp;
double delunm[5];
tmp = 1. / (omega * (2. - omega));
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP13896>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
}
rhs();
l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsdnm);
timer_clear(1);
timer_start(1);
for(istep = 1;
istep <= itmax;
istep++){
/*CompoundStmt:13989*/
if (istep % 20 == 0 || istep == itmax || istep == 1) {
    //#pragma omp master
        printf(" Time step %4d\n", istep);
}
if (affinMaskRes)
{
myDARTSRuntime->run(launch<TP13995>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (double*)&((tmp))));
}
if (istep % inorm == 0) {
    l2norm(nx0, ny0, nz0, ist, iend, jst, jend, delunm);
}
/*OMPMasterDirective:14144*/
printf("l2norm(1) complete...\n");
rhs();
/*OMPMasterDirective:14147*/
printf("rhs complete...\n");
if ((istep % inorm == 0) || (istep == itmax)) {
    l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsdnm);
}
/*OMPMasterDirective:14153*/
{
/*CompoundStmt:14154*/
printf("l2norm(2) complete...\n");
printf("%.2f %.2f %.2f %.2f %.2f\n", rsdnm[0], rsdnm[1], rsdnm[2], rsdnm[3], rsdnm[4]);
}
if ((rsdnm[0] < tolrsd[0]) && (rsdnm[1] < tolrsd[1]) && (rsdnm[2] < tolrsd[2]) && (rsdnm[3] < tolrsd[3]) && (rsdnm[4] < tolrsd[4])) {
    exit(1);
}
}
timer_stop(1);
maxtime = timer_read(1);
}
/*Function: verify, ID: 17*/
static void verify(double xcr[5], double xce[5], double xci, char *class_is, boolean *verified) {
/*verify:17*/
/*CompoundStmt:14169*/
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
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP192::_checkInCodelets194::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->nthreads_darts194 = (this->inputsTPParent->nthreads_darts192)/*OMP_SHARED - VAR INLINED*/;

/*printing node 195: BinaryOperator*/
(*(this->inputsTPParent->nthreads_darts194)) = omp_get_num_threads();
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/
myTP->controlTPParent->TPParent->barrierCodelets192[0].decDep();
}
else
{
/*Find and signal the next codelet*/
myTP->TPParent->barrierCodelets192[0].decDep();
}
}
TP192::TP192(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int *in_nthreads):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),nthreads_darts192(in_nthreads)/*OMP_SHARED - INPUT*/, TP194_alreadyLaunched(0) ,barrierCodelets192(new _barrierCodelets192[1]) ,checkInCodelets194(new _checkInCodelets194[this->numThreads]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets192[0] = _barrierCodelets192(ompNumThreads,ompNumThreads,this, 0);
_checkInCodelets194 * checkInCodelets194Ptr = (this->checkInCodelets194);
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets194Ptr) = _checkInCodelets194(1,1,this,codeletCounter);
(*checkInCodelets194Ptr).decDep();
checkInCodelets194Ptr++;
}
}
TP192::~TP192(){
delete [] barrierCodelets192;
delete [] checkInCodelets194;
}
/*TP1: TP_blts*/
void TP1::_checkInCodelets248::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*region 248 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP248;
if(idx < myTP->TPsToUse248){
if (!__sync_val_compare_and_swap (&(myTP->TP248_alreadyLaunched[idx]), 0, 1)){
int range = abs ((this->inputsTPParent->iend_darts1[this->getID()]) - (this->inputsTPParent->ist_darts1[this->getID()])) / 1;
int rangePerCodelet = range / myTP->TPsToUse248;
int minIteration = min<int >((this->inputsTPParent->iend_darts1[this->getID()]), (this->inputsTPParent->ist_darts1[this->getID()]));
int remainderRange = range % myTP->TPsToUse248;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if((this->inputsTPParent->ist_darts1[this->getID()]) < (this->inputsTPParent->iend_darts1[this->getID()]))
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse248 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse248 - 1)
{
lastIteration = (this->inputsTPParent->iend_darts1[this->getID()]);
}
#if USEINVOKE == 1
invoke < TP248 > (myTP, myTP->codeletsPerTP248 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP248Ptr[idx]));
#else
place < TP248 > (idx, myTP, myTP->codeletsPerTP248 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP248Ptr[idx]));
#endif
}else{
if (myTP->TP248Ptr[idx] != nullptr){
myTP->TP248Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP1::_barrierCodelets248::fire(void)
{
TP1* myTP =  static_cast<TP1*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets352[codeletsCounter].decDep();
}
}
}
void TP1::_checkInCodelets352::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/

/*printing node 353: CallExpr*/
printf("part1\n");
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 352 nextRegion: 354 */
myTP->controlTPParent->checkInCodelets354[this->getID()].decDep();
}
else
{
/*Signaling next codelet region: 352 nextRegion: 354 */
myTP->checkInCodelets354[this->getID()].decDep();
}
}
void TP1::_checkInCodelets354::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->iend_darts354= &(this->inputsTPParent->iend_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->ist_darts354= &(this->inputsTPParent->ist_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->jend_darts354= &(this->inputsTPParent->jend_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->jst_darts354= &(this->inputsTPParent->jst_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->k_darts354= &(this->inputsTPParent->k_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->omega_darts354= &(this->inputsTPParent->omega_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->tmp_darts354= &(this->inputsTPParent->tmp_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->tmp1_darts354= &(this->inputsTPParent->tmp1_darts1[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;

/*printing node 355: ForStmt*/
{
/*Loop's init*/
(this->inputsTPParent->i_darts354) = (*(this->inputsTPParent->ist_darts354));
int i_darts_counter_temp354 = (this->inputsTPParent->i_darts354);
for(;i_darts_counter_temp354 <= (*(this->inputsTPParent->iend_darts354));i_darts_counter_temp354++){
this->inputsTPParent->id_darts354 = omp_get_thread_num();
printf("Thread #%d: loopi %d start\n", (this->inputsTPParent->id_darts354), i_darts_counter_temp354);
if(i_darts_counter_temp354 != (*(this->inputsTPParent->ist_darts354)))
{
while(flag[i_darts_counter_temp354 - 1] == 0){
}
}
if(i_darts_counter_temp354 != (*(this->inputsTPParent->iend_darts354)))
{
while(flag[i_darts_counter_temp354] == 1){
}
}
{
/*Loop's init*/
(this->inputsTPParent->j_darts354) = (*(this->inputsTPParent->jst_darts354));
int j_darts_counter_temp354 = (this->inputsTPParent->j_darts354);
for(;j_darts_counter_temp354 <= (*(this->inputsTPParent->jend_darts354));j_darts_counter_temp354++){
{
/*Loop's init*/
(this->inputsTPParent->m_darts354) = 0;
int m_darts_counter_temp354 = (this->inputsTPParent->m_darts354);
for(;m_darts_counter_temp354 < 5;m_darts_counter_temp354++){
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][m_darts_counter_temp354] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][m_darts_counter_temp354] - (*(this->inputsTPParent->omega_darts354)) * (b[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][0] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354][j_darts_counter_temp354 - 1][0] + c[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][0] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354 - 1][j_darts_counter_temp354][0] + b[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][1] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354][j_darts_counter_temp354 - 1][1] + c[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][1] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354 - 1][j_darts_counter_temp354][1] + b[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][2] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354][j_darts_counter_temp354 - 1][2] + c[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][2] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354 - 1][j_darts_counter_temp354][2] + b[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][3] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354][j_darts_counter_temp354 - 1][3] + c[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][3] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354 - 1][j_darts_counter_temp354][3] + b[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][4] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354][j_darts_counter_temp354 - 1][4] + c[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][4] * rsd[(*(this->inputsTPParent->k_darts354))][i_darts_counter_temp354 - 1][j_darts_counter_temp354][4]);
}
(this->inputsTPParent->m_darts354) = m_darts_counter_temp354;
}
{
/*Loop's init*/
(this->inputsTPParent->m_darts354) = 0;
int m_darts_counter_temp354 = (this->inputsTPParent->m_darts354);
for(;m_darts_counter_temp354 < 5;m_darts_counter_temp354++){
tmat[m_darts_counter_temp354][0] = d[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][0];
tmat[m_darts_counter_temp354][1] = d[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][1];
tmat[m_darts_counter_temp354][2] = d[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][2];
tmat[m_darts_counter_temp354][3] = d[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][3];
tmat[m_darts_counter_temp354][4] = d[i_darts_counter_temp354][j_darts_counter_temp354][m_darts_counter_temp354][4];
}
(this->inputsTPParent->m_darts354) = m_darts_counter_temp354;
}
(*(this->inputsTPParent->tmp1_darts354)) = 1. / tmat[0][0];
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[1][0];
tmat[1][1] = tmat[1][1] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][1];
tmat[1][2] = tmat[1][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][2];
tmat[1][3] = tmat[1][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][3];
tmat[1][4] = tmat[1][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[2][0];
tmat[2][1] = tmat[2][1] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][1];
tmat[2][2] = tmat[2][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][2];
tmat[2][3] = tmat[2][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][3];
tmat[2][4] = tmat[2][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[3][0];
tmat[3][1] = tmat[3][1] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][1];
tmat[3][2] = tmat[3][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[4][0];
tmat[4][1] = tmat[4][1] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][1];
tmat[4][2] = tmat[4][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[0][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp1_darts354)) = 1. / tmat[1][1];
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[2][1];
tmat[2][2] = tmat[2][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][2];
tmat[2][3] = tmat[2][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][3];
tmat[2][4] = tmat[2][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[3][1];
tmat[3][2] = tmat[3][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[4][1];
tmat[4][2] = tmat[4][2] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[1][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp1_darts354)) = 1. / tmat[2][2];
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[3][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[2][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[2][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[4][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts354)) * tmat[2][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[2][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] * (*(this->inputsTPParent->tmp_darts354));
(*(this->inputsTPParent->tmp1_darts354)) = 1. / tmat[3][3];
(*(this->inputsTPParent->tmp_darts354)) = (*(this->inputsTPParent->tmp1_darts354)) * tmat[4][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts354)) * tmat[3][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] - rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] * (*(this->inputsTPParent->tmp_darts354));
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4] / tmat[4][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - tmat[3][4] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] / tmat[3][3];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] - tmat[2][3] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - tmat[2][4] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] / tmat[2][2];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] - tmat[1][2] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] - tmat[1][3] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - tmat[1][4] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] / tmat[1][1];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] - tmat[0][1] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][1] - tmat[0][2] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][2] - tmat[0][3] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][3] - tmat[0][4] * rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][4];
rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] = rsd[i_darts_counter_temp354][j_darts_counter_temp354][(*(this->inputsTPParent->k_darts354))][0] / tmat[0][0];
}
(this->inputsTPParent->j_darts354) = j_darts_counter_temp354;
}
if(i_darts_counter_temp354 != (*(this->inputsTPParent->ist_darts354)))
{
flag[i_darts_counter_temp354 - 1] = 0;
}
if(i_darts_counter_temp354 != (*(this->inputsTPParent->iend_darts354)))
{
flag[i_darts_counter_temp354] = 1;
}
}
(this->inputsTPParent->i_darts354) = i_darts_counter_temp354;
}
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsblts[this->getID()]->decDep();
}
else
{
/*Find and signal the next codelet*/

myTP->nextCodeletsblts[this->getID()]->decDep();
}
}
TP1::TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0):ompTP(in_numThreads, in_mainCodeletID), ptrToThisFunctionTP(in_ptrToThisFunctionTP), inputsTPParent(this), controlTPParent(this), nextCodeletsblts( new Codelet*[in_numThreads])
, nextSyncCodeletsblts( new Codelet*[in_numThreads])
,nx_darts1(new int[this->numThreads])
,ny_darts1(new int[this->numThreads])
,nz_darts1(new int[this->numThreads])
,k_darts1(new int[this->numThreads])
,omega_darts1(new double[this->numThreads])
,ist_darts1(new int[this->numThreads])
,iend_darts1(new int[this->numThreads])
,jst_darts1(new int[this->numThreads])
,jend_darts1(new int[this->numThreads])
,nx0_darts1(new int[this->numThreads])
,ny0_darts1(new int[this->numThreads])
,i_darts1(new int[this->numThreads])
,j_darts1(new int[this->numThreads])
,m_darts1(new int[this->numThreads])
,tmp_darts1(new double[this->numThreads])
,tmp1_darts1(new double[this->numThreads])
, TP248Ptr(new TP248 *[NUMTPS248]), TP248_alreadyLaunched(new size_t [NUMTPS248]), numTPsSet248(0), numTPsReady248(0), TPsToUse248(NUMTPS248), codeletsPerTP248(this->numThreads/NUMTPS248), totalCodelets248(this->TPsToUse248*this->codeletsPerTP248), TP352_alreadyLaunched(0), TP354_alreadyLaunched(0) ,checkInCodelets248(new _checkInCodelets248[this->numThreads]) ,barrierCodelets248(new _barrierCodelets248[1]) ,checkInCodelets352(new _checkInCodelets352[this->numThreads]) ,checkInCodelets354(new _checkInCodelets354[this->numThreads]){
barrierCodelets248[0] = _barrierCodelets248(NUMTPS248,NUMTPS248,this, 0);
_checkInCodelets354 * checkInCodelets354Ptr = (this->checkInCodelets354);
_checkInCodelets352 * checkInCodelets352Ptr = (this->checkInCodelets352);
_checkInCodelets248 * checkInCodelets248Ptr = (this->checkInCodelets248);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets248);
#endif
for(int i=0; i<NUMTPS248; i++)
{
TP248Ptr[i] = nullptr;
TP248_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++)
{
(*checkInCodelets354Ptr) = _checkInCodelets354(1,1,this,codeletCounter);
checkInCodelets354Ptr++;
(*checkInCodelets352Ptr) = _checkInCodelets352(1,1,this,codeletCounter);
checkInCodelets352Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets248Ptr) = _checkInCodelets248(2,1,this,codeletCounter);
#else
(*checkInCodelets248Ptr) = _checkInCodelets248(1,1,this,codeletCounter);
#endif
(*checkInCodelets248Ptr).decDep();
checkInCodelets248Ptr++;
}
if(this->numThreads == 1){
this->nextCodeletsblts[0] = in_mainNextCodelet;
this->nextSyncCodeletsblts[0] = in_mainSyncCodelet;
this->nx_darts1[0]= in_nx;
this->ny_darts1[0]= in_ny;
this->nz_darts1[0]= in_nz;
this->k_darts1[0]= in_k;
this->omega_darts1[0]= in_omega;
this->ist_darts1[0]= in_ist;
this->iend_darts1[0]= in_iend;
this->jst_darts1[0]= in_jst;
this->jend_darts1[0]= in_jend;
this->nx0_darts1[0]= in_nx0;
this->ny0_darts1[0]= in_ny0;
this->availableCodelets[0] = 1;
}
else
{
this->nx_darts1[this->mainCodeletID]= in_nx;
this->ny_darts1[this->mainCodeletID]= in_ny;
this->nz_darts1[this->mainCodeletID]= in_nz;
this->k_darts1[this->mainCodeletID]= in_k;
this->omega_darts1[this->mainCodeletID]= in_omega;
this->ist_darts1[this->mainCodeletID]= in_ist;
this->iend_darts1[this->mainCodeletID]= in_iend;
this->jst_darts1[this->mainCodeletID]= in_jst;
this->jend_darts1[this->mainCodeletID]= in_jend;
this->nx0_darts1[this->mainCodeletID]= in_nx0;
this->ny0_darts1[this->mainCodeletID]= in_ny0;
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
TP1::~TP1(){
delete [] checkInCodelets354;
delete [] checkInCodelets352;
delete [] barrierCodelets248;
delete [] checkInCodelets248;
delete [] nextCodeletsblts;
delete [] nextSyncCodeletsblts;
delete [] nx_darts1;
delete [] ny_darts1;
delete [] nz_darts1;
delete [] k_darts1;
delete [] omega_darts1;
delete [] ist_darts1;
delete [] iend_darts1;
delete [] jst_darts1;
delete [] jend_darts1;
delete [] nx0_darts1;
delete [] ny0_darts1;
delete [] i_darts1;
delete [] j_darts1;
delete [] m_darts1;
delete [] tmp_darts1;
delete [] tmp1_darts1;
}
void TP1::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID){
this->nx_darts1[codeletID]= in_nx;
this->ny_darts1[codeletID]= in_ny;
this->nz_darts1[codeletID]= in_nz;
this->k_darts1[codeletID]= in_k;
this->omega_darts1[codeletID]= in_omega;
this->ist_darts1[codeletID]= in_ist;
this->iend_darts1[codeletID]= in_iend;
this->jst_darts1[codeletID]= in_jst;
this->jend_darts1[codeletID]= in_jend;
this->nx0_darts1[codeletID]= in_nx0;
this->ny0_darts1[codeletID]= in_ny0;
}
/*TP248: OMPForDirective*/
void TP248::_barrierCodelets248::fire(void)
{
TP248* myTP = static_cast<TP248*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets248[0].decDep ();
}
bool TP248::requestNewRangeIterations248(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 1*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet248 * codeletID;
int tempEndRange   = rangePerCodelet248 * (codeletID + 1);
if (remainderRange248 != 0)
{
if (codeletID < (uint32_t)remainderRange248)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange248;
tempEndRange += remainderRange248;
}
}
tempStartRange = tempStartRange*1 + minIteration248;
tempEndRange = tempEndRange*1 + minIteration248;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration248 < lastIteration248)
{
(this->inputsTPParent->i_darts248[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts248[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration248;
}
}
return isThereNewIteration;
}
void TP248::_checkInCodelets249::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iend_darts248[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iend_darts1[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->ist_darts248[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->ist_darts1[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jend_darts248[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jend_darts1[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jst_darts248[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jst_darts1[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->k_darts248[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->k_darts1[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->omega_darts248[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->omega_darts1[this->getID()]);

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
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts248[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jend = &(this->inputsTPParent->jend_darts248[this->getLocalID()]);
(void)jend/*OMP_SHARED_PRIVATE*/;
int** jst = &(this->inputsTPParent->jst_darts248[this->getLocalID()]);
(void)jst/*OMP_SHARED_PRIVATE*/;
int** k = &(this->inputsTPParent->k_darts248[this->getLocalID()]);
(void)k/*OMP_SHARED_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts248[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** omega = &(this->inputsTPParent->omega_darts248[this->getLocalID()]);
(void)omega/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations248((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets248[0].decDep();
return;
}
for (int i_darts_counter_temp248 = (*i);i_darts_counter_temp248<=endRange && i_darts_counter_temp248<=this->inputsTPParent->lastIteration248;i_darts_counter_temp248++)
{
{
{
/*Loop's init*/
(*j) = (*(*jst));
int j_darts_counter_temp248 = (*j);
for(;j_darts_counter_temp248 <= (*(*jend));j_darts_counter_temp248++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp248 = (*m);
for(;m_darts_counter_temp248 < 5;m_darts_counter_temp248++){
rsd[(i_darts_counter_temp248)][j_darts_counter_temp248][(*(*k))][m_darts_counter_temp248] = rsd[(i_darts_counter_temp248)][j_darts_counter_temp248][(*(*k))][m_darts_counter_temp248] - (*(*omega)) * (a[(i_darts_counter_temp248)][j_darts_counter_temp248][m_darts_counter_temp248][0] * rsd[(*(*k)) - 1][(i_darts_counter_temp248)][j_darts_counter_temp248][0] + a[(i_darts_counter_temp248)][j_darts_counter_temp248][m_darts_counter_temp248][1] * rsd[(*(*k)) - 1][(i_darts_counter_temp248)][j_darts_counter_temp248][1] + a[(i_darts_counter_temp248)][j_darts_counter_temp248][m_darts_counter_temp248][2] * rsd[(*(*k)) - 1][(i_darts_counter_temp248)][j_darts_counter_temp248][2] + a[(i_darts_counter_temp248)][j_darts_counter_temp248][m_darts_counter_temp248][3] * rsd[(*(*k)) - 1][(i_darts_counter_temp248)][j_darts_counter_temp248][3] + a[(i_darts_counter_temp248)][j_darts_counter_temp248][m_darts_counter_temp248][4] * rsd[(*(*k)) - 1][(i_darts_counter_temp248)][j_darts_counter_temp248][4]);
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
TP248::TP248(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration, int in_lastIteration, TP248** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts248(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iend_darts248(new int*[this->numThreads]),ist_darts248(new int*[this->numThreads]),j_darts248(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jend_darts248(new int*[this->numThreads]),jst_darts248(new int*[this->numThreads]),k_darts248(new int*[this->numThreads]),m_darts248(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,omega_darts248(new double*[this->numThreads]), initIteration248(in_initIteration), lastIteration248(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets248(new _barrierCodelets248[1]) ,checkInCodelets249(new _checkInCodelets249[this->numThreads]){
/*Initialize the loop parameters*/
range248 = abs (lastIteration248 - initIteration248) / 1;
rangePerCodelet248 = range248 / numThreads;
minIteration248 = min<int>(lastIteration248, initIteration248);
remainderRange248 = range248 % numThreads;
/*Initialize inputs and vars.*/
this->iend_darts248 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->ist_darts248 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jend_darts248 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jst_darts248 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->k_darts248 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->omega_darts248 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets248[0] = _barrierCodelets248(this->numThreads,this->numThreads,this, 0);
_checkInCodelets249 * checkInCodelets249Ptr = (this->checkInCodelets249);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets249);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets249Ptr) = _checkInCodelets249(2,1,this,codeletCounter);
#else
(*checkInCodelets249Ptr) = _checkInCodelets249(1,1,this,codeletCounter);
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
this->checkInCodelets249[localID].setID (codeletID);
this->checkInCodelets249[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets249[localID + this->baseNumThreads * i] = _checkInCodelets249(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets249[localID + this->baseNumThreads * i] = _checkInCodelets249(1, 1, this, localID + this->baseNumThreads * i);
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
TP248::~TP248(){
delete [] iend_darts248;
delete [] ist_darts248;
delete [] jend_darts248;
delete [] jst_darts248;
delete [] k_darts248;
delete [] omega_darts248;
delete [] barrierCodelets248;
delete [] checkInCodelets249;
}
/*TP2: TP_buts*/
void TP2::_checkInCodelets1218::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*region 1218 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP1218;
if(idx < myTP->TPsToUse1218){
if (!__sync_val_compare_and_swap (&(myTP->TP1218_alreadyLaunched[idx]), 0, 1)){
int range = abs ((this->inputsTPParent->ist_darts2[this->getID()]) - (this->inputsTPParent->iend_darts2[this->getID()])) / 1;
int rangePerCodelet = range / myTP->TPsToUse1218;
int minIteration = min<int >((this->inputsTPParent->ist_darts2[this->getID()]), (this->inputsTPParent->iend_darts2[this->getID()]));
int remainderRange = range % myTP->TPsToUse1218;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if((this->inputsTPParent->iend_darts2[this->getID()]) < (this->inputsTPParent->ist_darts2[this->getID()]))
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == 0)
{
lastIteration = lastIteration - 1;
}
if(idx == myTP->TPsToUse1218 - 1)
{
lastIteration = (this->inputsTPParent->ist_darts2[this->getID()]);
}
#if USEINVOKE == 1
invoke < TP1218 > (myTP, myTP->codeletsPerTP1218 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP1218Ptr[idx]));
#else
place < TP1218 > (idx, myTP, myTP->codeletsPerTP1218 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP1218Ptr[idx]));
#endif
}else{
if (myTP->TP1218Ptr[idx] != nullptr){
myTP->TP1218Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
else
{
/*Signaling next codelet region: 1218 nextRegion: 1320 */
myTP->controlTPParent->checkInCodelets1320[this->getID()].decDep();
}
}
void TP2::_checkInCodelets1320::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->iend_darts1320= &(this->inputsTPParent->iend_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->ist_darts1320= &(this->inputsTPParent->ist_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->jend_darts1320= &(this->inputsTPParent->jend_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->jst_darts1320= &(this->inputsTPParent->jst_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->k_darts1320= &(this->inputsTPParent->k_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->omega_darts1320= &(this->inputsTPParent->omega_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->tmp_darts1320= &(this->inputsTPParent->tmp_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;
this->inputsTPParent->tmp1_darts1320= &(this->inputsTPParent->tmp1_darts2[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;

/*printing node 1321: ForStmt*/
{
/*Loop's init*/
(this->inputsTPParent->i_darts1320) = (*(this->inputsTPParent->iend_darts1320));
int i_darts_counter_temp1320 = (this->inputsTPParent->i_darts1320);
for(;i_darts_counter_temp1320 >= (*(this->inputsTPParent->ist_darts1320));i_darts_counter_temp1320--){
if(i_darts_counter_temp1320 != (*(this->inputsTPParent->iend_darts1320)))
{
while(flag[i_darts_counter_temp1320 + 1] == 0){
}
}
if(i_darts_counter_temp1320 != (*(this->inputsTPParent->ist_darts1320)))
{
while(flag[i_darts_counter_temp1320] == 1){
}
}
{
/*Loop's init*/
(this->inputsTPParent->j_darts1320) = (*(this->inputsTPParent->jend_darts1320));
int j_darts_counter_temp1320 = (this->inputsTPParent->j_darts1320);
for(;j_darts_counter_temp1320 >= (*(this->inputsTPParent->jst_darts1320));j_darts_counter_temp1320--){
{
/*Loop's init*/
(this->inputsTPParent->m_darts1320) = 0;
int m_darts_counter_temp1320 = (this->inputsTPParent->m_darts1320);
for(;m_darts_counter_temp1320 < 5;m_darts_counter_temp1320++){
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320] + (*(this->inputsTPParent->omega_darts1320)) * (b[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][0] * rsd[i_darts_counter_temp1320][j_darts_counter_temp1320 + 1][(*(this->inputsTPParent->k_darts1320))][0] + a[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][0] * rsd[i_darts_counter_temp1320 + 1][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][0] + b[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][1] * rsd[i_darts_counter_temp1320][j_darts_counter_temp1320 + 1][(*(this->inputsTPParent->k_darts1320))][1] + a[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][1] * rsd[i_darts_counter_temp1320 + 1][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][1] + b[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][2] * rsd[i_darts_counter_temp1320][j_darts_counter_temp1320 + 1][(*(this->inputsTPParent->k_darts1320))][2] + a[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][2] * rsd[i_darts_counter_temp1320 + 1][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][2] + b[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][3] * rsd[i_darts_counter_temp1320][j_darts_counter_temp1320 + 1][(*(this->inputsTPParent->k_darts1320))][3] + a[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][3] * rsd[i_darts_counter_temp1320 + 1][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][3] + b[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][4] * rsd[i_darts_counter_temp1320][j_darts_counter_temp1320 + 1][(*(this->inputsTPParent->k_darts1320))][4] + a[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][4] * rsd[i_darts_counter_temp1320 + 1][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][4]);
}
(this->inputsTPParent->m_darts1320) = m_darts_counter_temp1320;
}
{
/*Loop's init*/
(this->inputsTPParent->m_darts1320) = 0;
int m_darts_counter_temp1320 = (this->inputsTPParent->m_darts1320);
for(;m_darts_counter_temp1320 < 5;m_darts_counter_temp1320++){
tmat[m_darts_counter_temp1320][0] = d[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][0];
tmat[m_darts_counter_temp1320][1] = d[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][1];
tmat[m_darts_counter_temp1320][2] = d[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][2];
tmat[m_darts_counter_temp1320][3] = d[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][3];
tmat[m_darts_counter_temp1320][4] = d[i_darts_counter_temp1320][j_darts_counter_temp1320][m_darts_counter_temp1320][4];
}
(this->inputsTPParent->m_darts1320) = m_darts_counter_temp1320;
}
(*(this->inputsTPParent->tmp1_darts1320)) = 1. / tmat[0][0];
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[1][0];
tmat[1][1] = tmat[1][1] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][1];
tmat[1][2] = tmat[1][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][2];
tmat[1][3] = tmat[1][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][3];
tmat[1][4] = tmat[1][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[2][0];
tmat[2][1] = tmat[2][1] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][1];
tmat[2][2] = tmat[2][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][2];
tmat[2][3] = tmat[2][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][3];
tmat[2][4] = tmat[2][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[3][0];
tmat[3][1] = tmat[3][1] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][1];
tmat[3][2] = tmat[3][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[4][0];
tmat[4][1] = tmat[4][1] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][1];
tmat[4][2] = tmat[4][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[0][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp1_darts1320)) = 1. / tmat[1][1];
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[2][1];
tmat[2][2] = tmat[2][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][2];
tmat[2][3] = tmat[2][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][3];
tmat[2][4] = tmat[2][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[3][1];
tmat[3][2] = tmat[3][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[4][1];
tmat[4][2] = tmat[4][2] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[1][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp1_darts1320)) = 1. / tmat[2][2];
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[3][2];
tmat[3][3] = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[2][3];
tmat[3][4] = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[2][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[4][2];
tmat[4][3] = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[2][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[2][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] * (*(this->inputsTPParent->tmp_darts1320));
(*(this->inputsTPParent->tmp1_darts1320)) = 1. / tmat[3][3];
(*(this->inputsTPParent->tmp_darts1320)) = (*(this->inputsTPParent->tmp1_darts1320)) * tmat[4][3];
tmat[4][4] = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1320)) * tmat[3][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] * (*(this->inputsTPParent->tmp_darts1320));
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4] / tmat[4][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tmat[3][4] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] / tmat[3][3];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] - tmat[2][3] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tmat[2][4] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] / tmat[2][2];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] - tmat[1][2] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] - tmat[1][3] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tmat[1][4] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] / tmat[1][1];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] - tmat[0][1] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1] - tmat[0][2] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2] - tmat[0][3] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3] - tmat[0][4] * tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4];
tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] = tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0] / tmat[0][0];
rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][0] = rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][0] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][0];
rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][1] = rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][1] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][1];
rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][2] = rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][2] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][2];
rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][3] = rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][3] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][3];
rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][4] = rsd[i_darts_counter_temp1320][j_darts_counter_temp1320][(*(this->inputsTPParent->k_darts1320))][4] - tv[i_darts_counter_temp1320][j_darts_counter_temp1320][4];
}
(this->inputsTPParent->j_darts1320) = j_darts_counter_temp1320;
}
if(i_darts_counter_temp1320 != (*(this->inputsTPParent->iend_darts1320)))
{
flag[i_darts_counter_temp1320 + 1] = 0;
}
if(i_darts_counter_temp1320 != (*(this->inputsTPParent->ist_darts1320)))
{
flag[i_darts_counter_temp1320] = 1;
}
}
(this->inputsTPParent->i_darts1320) = i_darts_counter_temp1320;
}
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsbuts[this->getID()]->decDep();
}
else
{
/*Find and signal the next codelet*/

myTP->nextCodeletsbuts[this->getID()]->decDep();
}
}
TP2::TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0):ompTP(in_numThreads, in_mainCodeletID), ptrToThisFunctionTP(in_ptrToThisFunctionTP), inputsTPParent(this), controlTPParent(this), nextCodeletsbuts( new Codelet*[in_numThreads])
, nextSyncCodeletsbuts( new Codelet*[in_numThreads])
,nx_darts2(new int[this->numThreads])
,ny_darts2(new int[this->numThreads])
,nz_darts2(new int[this->numThreads])
,k_darts2(new int[this->numThreads])
,omega_darts2(new double[this->numThreads])
,ist_darts2(new int[this->numThreads])
,iend_darts2(new int[this->numThreads])
,jst_darts2(new int[this->numThreads])
,jend_darts2(new int[this->numThreads])
,nx0_darts2(new int[this->numThreads])
,ny0_darts2(new int[this->numThreads])
,i_darts2(new int[this->numThreads])
,j_darts2(new int[this->numThreads])
,m_darts2(new int[this->numThreads])
,tmp_darts2(new double[this->numThreads])
,tmp1_darts2(new double[this->numThreads])
, TP1218Ptr(new TP1218 *[NUMTPS1218]), TP1218_alreadyLaunched(new size_t [NUMTPS1218]), numTPsSet1218(0), numTPsReady1218(0), TPsToUse1218(NUMTPS1218), codeletsPerTP1218(this->numThreads/NUMTPS1218), totalCodelets1218(this->TPsToUse1218*this->codeletsPerTP1218), TP1320_alreadyLaunched(0) ,checkInCodelets1218(new _checkInCodelets1218[this->numThreads]) ,checkInCodelets1320(new _checkInCodelets1320[this->numThreads]){
_checkInCodelets1320 * checkInCodelets1320Ptr = (this->checkInCodelets1320);
_checkInCodelets1218 * checkInCodelets1218Ptr = (this->checkInCodelets1218);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets1218);
#endif
for(int i=0; i<NUMTPS1218; i++)
{
TP1218Ptr[i] = nullptr;
TP1218_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++)
{
(*checkInCodelets1320Ptr) = _checkInCodelets1320(1,1,this,codeletCounter);
checkInCodelets1320Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets1218Ptr) = _checkInCodelets1218(2,1,this,codeletCounter);
#else
(*checkInCodelets1218Ptr) = _checkInCodelets1218(1,1,this,codeletCounter);
#endif
(*checkInCodelets1218Ptr).decDep();
checkInCodelets1218Ptr++;
}
if(this->numThreads == 1){
this->nextCodeletsbuts[0] = in_mainNextCodelet;
this->nextSyncCodeletsbuts[0] = in_mainSyncCodelet;
this->nx_darts2[0]= in_nx;
this->ny_darts2[0]= in_ny;
this->nz_darts2[0]= in_nz;
this->k_darts2[0]= in_k;
this->omega_darts2[0]= in_omega;
this->ist_darts2[0]= in_ist;
this->iend_darts2[0]= in_iend;
this->jst_darts2[0]= in_jst;
this->jend_darts2[0]= in_jend;
this->nx0_darts2[0]= in_nx0;
this->ny0_darts2[0]= in_ny0;
this->availableCodelets[0] = 1;
}
else
{
this->nx_darts2[this->mainCodeletID]= in_nx;
this->ny_darts2[this->mainCodeletID]= in_ny;
this->nz_darts2[this->mainCodeletID]= in_nz;
this->k_darts2[this->mainCodeletID]= in_k;
this->omega_darts2[this->mainCodeletID]= in_omega;
this->ist_darts2[this->mainCodeletID]= in_ist;
this->iend_darts2[this->mainCodeletID]= in_iend;
this->jst_darts2[this->mainCodeletID]= in_jst;
this->jend_darts2[this->mainCodeletID]= in_jend;
this->nx0_darts2[this->mainCodeletID]= in_nx0;
this->ny0_darts2[this->mainCodeletID]= in_ny0;
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
TP2::~TP2(){
delete [] checkInCodelets1320;
delete [] checkInCodelets1218;
delete [] nextCodeletsbuts;
delete [] nextSyncCodeletsbuts;
delete [] nx_darts2;
delete [] ny_darts2;
delete [] nz_darts2;
delete [] k_darts2;
delete [] omega_darts2;
delete [] ist_darts2;
delete [] iend_darts2;
delete [] jst_darts2;
delete [] jend_darts2;
delete [] nx0_darts2;
delete [] ny0_darts2;
delete [] i_darts2;
delete [] j_darts2;
delete [] m_darts2;
delete [] tmp_darts2;
delete [] tmp1_darts2;
}
void TP2::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID){
this->nx_darts2[codeletID]= in_nx;
this->ny_darts2[codeletID]= in_ny;
this->nz_darts2[codeletID]= in_nz;
this->k_darts2[codeletID]= in_k;
this->omega_darts2[codeletID]= in_omega;
this->ist_darts2[codeletID]= in_ist;
this->iend_darts2[codeletID]= in_iend;
this->jst_darts2[codeletID]= in_jst;
this->jend_darts2[codeletID]= in_jend;
this->nx0_darts2[codeletID]= in_nx0;
this->ny0_darts2[codeletID]= in_ny0;
}
/*TP1218: OMPForDirective*/
bool TP1218::requestNewRangeIterations1218(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 1*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet1218 * codeletID;
int tempEndRange   = rangePerCodelet1218 * (codeletID + 1);
if (remainderRange1218 != 0)
{
if (codeletID < (uint32_t)remainderRange1218)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange1218;
tempEndRange += remainderRange1218;
}
}
tempStartRange = tempStartRange*1 + minIteration1218;
tempEndRange = tempEndRange*1 + minIteration1218;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration1218 < lastIteration1218)
{
(this->inputsTPParent->i_darts1218[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts1218[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == 0)
{
*endRange = *endRange - 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration1218;
}
}
return isThereNewIteration;
}
void TP1218::_checkInCodelets1219::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iend_darts1218[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iend_darts2[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->ist_darts1218[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->ist_darts2[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jend_darts1218[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jend_darts2[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jst_darts1218[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jst_darts2[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->k_darts1218[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->k_darts2[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->omega_darts1218[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->omega_darts2[this->getID()]);

/*printing node 1219: ForStmt*/
/*var: i*/
/*var: iend*/
/*var: ist*/
/*var: j*/
/*var: jend*/
/*var: jst*/
/*var: k*/
/*var: m*/
/*var: omega*/
int* i = &(this->inputsTPParent->i_darts1218[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts1218[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jend = &(this->inputsTPParent->jend_darts1218[this->getLocalID()]);
(void)jend/*OMP_SHARED_PRIVATE*/;
int** jst = &(this->inputsTPParent->jst_darts1218[this->getLocalID()]);
(void)jst/*OMP_SHARED_PRIVATE*/;
int** k = &(this->inputsTPParent->k_darts1218[this->getLocalID()]);
(void)k/*OMP_SHARED_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts1218[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** omega = &(this->inputsTPParent->omega_darts1218[this->getLocalID()]);
(void)omega/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations1218((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Find and signal the next codelet*/
myTP->controlTPParent->TPParent->checkInCodelets1320[this->getID()].decDep();
return;
}
for (int i_darts_counter_temp1218 = (*i);i_darts_counter_temp1218>=endRange && i_darts_counter_temp1218>=this->inputsTPParent->lastIteration1218;i_darts_counter_temp1218--)
{
{
{
/*Loop's init*/
(*j) = (*(*jend));
int j_darts_counter_temp1218 = (*j);
for(;j_darts_counter_temp1218 >= (*(*jst));j_darts_counter_temp1218--){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp1218 = (*m);
for(;m_darts_counter_temp1218 < 5;m_darts_counter_temp1218++){
tv[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218] = (*(*omega)) * (c[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218][0] * rsd[(i_darts_counter_temp1218)][j_darts_counter_temp1218][(*(*k)) + 1][0] + c[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218][1] * rsd[(i_darts_counter_temp1218)][j_darts_counter_temp1218][(*(*k)) + 1][1] + c[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218][2] * rsd[(i_darts_counter_temp1218)][j_darts_counter_temp1218][(*(*k)) + 1][2] + c[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218][3] * rsd[(i_darts_counter_temp1218)][j_darts_counter_temp1218][(*(*k)) + 1][3] + c[(i_darts_counter_temp1218)][j_darts_counter_temp1218][m_darts_counter_temp1218][4] * rsd[(i_darts_counter_temp1218)][j_darts_counter_temp1218][(*(*k)) + 1][4]);
}
(*m) = m_darts_counter_temp1218;
}
}
(*j) = j_darts_counter_temp1218;
}
}
}
/*If this omp for has no barrier, 
check if all the codelets 
replicated from the same 
global ID has finished and 
signal the next codelet. 
Otherwise, return.*/
uint32_t completedMultCodelet = __sync_fetch_and_add(&(myTP->signalNextReady[this->getLocalID() % myTP->baseNumThreads]), 1);
if(completedMultCodelet < (uint32_t)(DARTS_CODELETS_MULT - 1))
return;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/
myTP->controlTPParent->TPParent->checkInCodelets1320[this->getID()].decDep();
}
TP1218::TP1218(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration, int in_lastIteration, TP1218** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts1218(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iend_darts1218(new int*[this->numThreads]),ist_darts1218(new int*[this->numThreads]),j_darts1218(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jend_darts1218(new int*[this->numThreads]),jst_darts1218(new int*[this->numThreads]),k_darts1218(new int*[this->numThreads]),m_darts1218(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,omega_darts1218(new double*[this->numThreads]), initIteration1218(in_initIteration), lastIteration1218(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT), signalNextReady(new int[baseNumThreads]) ,checkInCodelets1219(new _checkInCodelets1219[this->numThreads]){
/*Initialize the loop parameters*/
range1218 = abs (lastIteration1218 - initIteration1218) / 1;
rangePerCodelet1218 = range1218 / numThreads;
minIteration1218 = min<int>(lastIteration1218, initIteration1218);
remainderRange1218 = range1218 % numThreads;
/*Initialize inputs and vars.*/
this->iend_darts1218 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->ist_darts1218 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jend_darts1218 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jst_darts1218 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->k_darts1218 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->omega_darts1218 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
_checkInCodelets1219 * checkInCodelets1219Ptr = (this->checkInCodelets1219);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets1219);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets1219Ptr) = _checkInCodelets1219(2,1,this,codeletCounter);
#else
(*checkInCodelets1219Ptr) = _checkInCodelets1219(1,1,this,codeletCounter);
#endif
(*checkInCodelets1219Ptr).decDep();
checkInCodelets1219Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP1218::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets1219[localID].setID (codeletID);
this->checkInCodelets1219[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets1219[localID + this->baseNumThreads * i] = _checkInCodelets1219(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets1219[localID + this->baseNumThreads * i] = _checkInCodelets1219(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets1219[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets1219[localID + this->baseNumThreads * i].decDep();
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
TP1218::~TP1218(){
delete [] iend_darts1218;
delete [] ist_darts1218;
delete [] jend_darts1218;
delete [] jst_darts1218;
delete [] k_darts1218;
delete [] omega_darts1218;
delete [] checkInCodelets1219;
}
/*TP2204: OMPParallelDirective*/
void TP2204::_barrierCodelets2204::fire(void)
{
TP2204* myTP = static_cast<TP2204*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP2204::_checkInCodelets2206::fire(void)
{
/*Init the vars for this region*/

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

/*printing node 2217: DeclStmt*/

/*printing node 2218: DeclStmt*/

/*printing node 2219: DeclStmt*/

/*printing node 2220: DeclStmt*/

/*printing node 2221: DeclStmt*/

/*printing node 2222: BinaryOperator*/
(this->inputsTPParent->dsspm_darts2204[this->getID()]) = dssp;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 2206 nextRegion: 2223 */
myTP->controlTPParent->checkInCodelets2223[this->getID()].decDep();
}
void TP2204::_checkInCodelets2223::fire(void)
{
/*region 2223 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP2223;
if(idx < myTP->TPsToUse2223){
if (!__sync_val_compare_and_swap (&(myTP->TP2223_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse2223;
int minIteration = min<int >(nx, 0);
int remainderRange = range % myTP->TPsToUse2223;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse2223 - 1)
{
lastIteration = nx;
}
#if USEINVOKE == 1
invoke < TP2223 > (myTP, myTP->codeletsPerTP2223 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2223Ptr[idx]));
#else
place < TP2223 > (idx, myTP, myTP->codeletsPerTP2223 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2223Ptr[idx]));
#endif
}else{
if (myTP->TP2223Ptr[idx] != nullptr){
myTP->TP2223Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets2223::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets2274[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets2274::fire(void)
{
/*region 2274 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP2274;
if(idx < myTP->TPsToUse2274){
if (!__sync_val_compare_and_swap (&(myTP->TP2274_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse2274;
int minIteration = min<int >(nx, 0);
int remainderRange = range % myTP->TPsToUse2274;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse2274 - 1)
{
lastIteration = nx;
}
#if USEINVOKE == 1
invoke < TP2274 > (myTP, myTP->codeletsPerTP2274 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2274Ptr[idx]));
#else
place < TP2274 > (idx, myTP, myTP->codeletsPerTP2274 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2274Ptr[idx]));
#endif
}else{
if (myTP->TP2274Ptr[idx] != nullptr){
myTP->TP2274Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets2274::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets2406[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets2406::fire(void)
{

/*printing node 2406: BinaryOperator*/
(this->inputsTPParent->L1_darts2204[this->getID()]) = 0;

/*printing node 2407: BinaryOperator*/
(this->inputsTPParent->L2_darts2204[this->getID()]) = nx - 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 2406 nextRegion: 2409 */
myTP->controlTPParent->checkInCodelets2409[this->getID()].decDep();
}
void TP2204::_checkInCodelets2409::fire(void)
{
/*region 2409 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP2409;
if(idx < myTP->TPsToUse2409){
if (!__sync_val_compare_and_swap (&(myTP->TP2409_alreadyLaunched[idx]), 0, 1)){
int range = abs ((this->inputsTPParent->L2_darts2204[this->getID()]) - (this->inputsTPParent->L1_darts2204[this->getID()])) / 1;
int rangePerCodelet = range / myTP->TPsToUse2409;
int minIteration = min<int >((this->inputsTPParent->L2_darts2204[this->getID()]), (this->inputsTPParent->L1_darts2204[this->getID()]));
int remainderRange = range % myTP->TPsToUse2409;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if((this->inputsTPParent->L1_darts2204[this->getID()]) < (this->inputsTPParent->L2_darts2204[this->getID()]))
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse2409 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse2409 - 1)
{
lastIteration = (this->inputsTPParent->L2_darts2204[this->getID()]);
}
#if USEINVOKE == 1
invoke < TP2409 > (myTP, myTP->codeletsPerTP2409 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2409Ptr[idx]));
#else
place < TP2409 > (idx, myTP, myTP->codeletsPerTP2409 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2409Ptr[idx]));
#endif
}else{
if (myTP->TP2409Ptr[idx] != nullptr){
myTP->TP2409Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets2409::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets2558[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets2558::fire(void)
{
/*region 2558 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP2558;
if(idx < myTP->TPsToUse2558){
if (!__sync_val_compare_and_swap (&(myTP->TP2558_alreadyLaunched[idx]), 0, 1)){
int range = abs (jend - jst) / 1;
int rangePerCodelet = range / myTP->TPsToUse2558;
int minIteration = min<int >(jend, jst);
int remainderRange = range % myTP->TPsToUse2558;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(jst < jend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse2558 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse2558 - 1)
{
lastIteration = jend;
}
#if USEINVOKE == 1
invoke < TP2558 > (myTP, myTP->codeletsPerTP2558 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2558Ptr[idx]));
#else
place < TP2558 > (idx, myTP, myTP->codeletsPerTP2558 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP2558Ptr[idx]));
#endif
}else{
if (myTP->TP2558Ptr[idx] != nullptr){
myTP->TP2558Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets2558::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets3193[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets3193::fire(void)
{

/*printing node 3193: BinaryOperator*/
(this->inputsTPParent->L1_darts2204[this->getID()]) = 0;

/*printing node 3194: BinaryOperator*/
(this->inputsTPParent->L2_darts2204[this->getID()]) = ny - 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 3193 nextRegion: 3196 */
myTP->controlTPParent->checkInCodelets3196[this->getID()].decDep();
}
void TP2204::_checkInCodelets3196::fire(void)
{
/*region 3196 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP3196;
if(idx < myTP->TPsToUse3196){
if (!__sync_val_compare_and_swap (&(myTP->TP3196_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse3196;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse3196;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse3196 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse3196 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP3196 > (myTP, myTP->codeletsPerTP3196 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3196Ptr[idx]));
#else
place < TP3196 > (idx, myTP, myTP->codeletsPerTP3196 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3196Ptr[idx]));
#endif
}else{
if (myTP->TP3196Ptr[idx] != nullptr){
myTP->TP3196Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets3196::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets3345[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets3345::fire(void)
{
/*region 3345 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP3345;
if(idx < myTP->TPsToUse3345){
if (!__sync_val_compare_and_swap (&(myTP->TP3345_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse3345;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse3345;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse3345 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse3345 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP3345 > (myTP, myTP->codeletsPerTP3345 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3345Ptr[idx]));
#else
place < TP3345 > (idx, myTP, myTP->codeletsPerTP3345 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3345Ptr[idx]));
#endif
}else{
if (myTP->TP3345Ptr[idx] != nullptr){
myTP->TP3345Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets3345::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets3980[codeletsCounter].decDep();
}
}
}
void TP2204::_checkInCodelets3980::fire(void)
{
/*region 3980 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP3980;
if(idx < myTP->TPsToUse3980){
if (!__sync_val_compare_and_swap (&(myTP->TP3980_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse3980;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse3980;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse3980 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse3980 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP3980 > (myTP, myTP->codeletsPerTP3980 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3980Ptr[idx]));
#else
place < TP3980 > (idx, myTP, myTP->codeletsPerTP3980 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP3980Ptr[idx]));
#endif
}else{
if (myTP->TP3980Ptr[idx] != nullptr){
myTP->TP3980Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP2204::_barrierCodelets3980::fire(void)
{
TP2204* myTP =  static_cast<TP2204*>(myTP_);
myTP->TPParent->barrierCodelets2204[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets2204[0]));
}
TP2204::TP2204(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),L1_darts2204(new int[this->numThreads])/*VARIABLE*/,L2_darts2204(new int[this->numThreads])/*VARIABLE*/,dsspm_darts2204(new double[this->numThreads])/*VARIABLE*/,eta_darts2204(new double[this->numThreads])/*VARIABLE*/,i_darts2204(new int[this->numThreads])/*VARIABLE*/,iend1_darts2204(new int[this->numThreads])/*VARIABLE*/,iglob_darts2204(new int[this->numThreads])/*VARIABLE*/,ist1_darts2204(new int[this->numThreads])/*VARIABLE*/,j_darts2204(new int[this->numThreads])/*VARIABLE*/,jend1_darts2204(new int[this->numThreads])/*VARIABLE*/,jglob_darts2204(new int[this->numThreads])/*VARIABLE*/,jst1_darts2204(new int[this->numThreads])/*VARIABLE*/,k_darts2204(new int[this->numThreads])/*VARIABLE*/,m_darts2204(new int[this->numThreads])/*VARIABLE*/,q_darts2204(new double[this->numThreads])/*VARIABLE*/,tmp_darts2204(new double[this->numThreads])/*VARIABLE*/,u21_darts2204(new double[this->numThreads])/*VARIABLE*/,u21i_darts2204(new double[this->numThreads])/*VARIABLE*/,u21im1_darts2204(new double[this->numThreads])/*VARIABLE*/,u21j_darts2204(new double[this->numThreads])/*VARIABLE*/,u21jm1_darts2204(new double[this->numThreads])/*VARIABLE*/,u21k_darts2204(new double[this->numThreads])/*VARIABLE*/,u21km1_darts2204(new double[this->numThreads])/*VARIABLE*/,u31_darts2204(new double[this->numThreads])/*VARIABLE*/,u31i_darts2204(new double[this->numThreads])/*VARIABLE*/,u31im1_darts2204(new double[this->numThreads])/*VARIABLE*/,u31j_darts2204(new double[this->numThreads])/*VARIABLE*/,u31jm1_darts2204(new double[this->numThreads])/*VARIABLE*/,u31k_darts2204(new double[this->numThreads])/*VARIABLE*/,u31km1_darts2204(new double[this->numThreads])/*VARIABLE*/,u41_darts2204(new double[this->numThreads])/*VARIABLE*/,u41i_darts2204(new double[this->numThreads])/*VARIABLE*/,u41im1_darts2204(new double[this->numThreads])/*VARIABLE*/,u41j_darts2204(new double[this->numThreads])/*VARIABLE*/,u41jm1_darts2204(new double[this->numThreads])/*VARIABLE*/,u41k_darts2204(new double[this->numThreads])/*VARIABLE*/,u41km1_darts2204(new double[this->numThreads])/*VARIABLE*/,u51i_darts2204(new double[this->numThreads])/*VARIABLE*/,u51im1_darts2204(new double[this->numThreads])/*VARIABLE*/,u51j_darts2204(new double[this->numThreads])/*VARIABLE*/,u51jm1_darts2204(new double[this->numThreads])/*VARIABLE*/,u51k_darts2204(new double[this->numThreads])/*VARIABLE*/,u51km1_darts2204(new double[this->numThreads])/*VARIABLE*/,xi_darts2204(new double[this->numThreads])/*VARIABLE*/,zeta_darts2204(new double[this->numThreads])/*VARIABLE*/, TP2223Ptr(new TP2223 *[NUMTPS2223]), TP2223_alreadyLaunched(new size_t [NUMTPS2223]), numTPsSet2223(0), numTPsReady2223(0), TPsToUse2223(NUMTPS2223), codeletsPerTP2223(this->numThreads/NUMTPS2223), totalCodelets2223(this->TPsToUse2223*this->codeletsPerTP2223), TP2274Ptr(new TP2274 *[NUMTPS2274]), TP2274_alreadyLaunched(new size_t [NUMTPS2274]), numTPsSet2274(0), numTPsReady2274(0), TPsToUse2274(NUMTPS2274), codeletsPerTP2274(this->numThreads/NUMTPS2274), totalCodelets2274(this->TPsToUse2274*this->codeletsPerTP2274), TP2409Ptr(new TP2409 *[NUMTPS2409]), TP2409_alreadyLaunched(new size_t [NUMTPS2409]), numTPsSet2409(0), numTPsReady2409(0), TPsToUse2409(NUMTPS2409), codeletsPerTP2409(this->numThreads/NUMTPS2409), totalCodelets2409(this->TPsToUse2409*this->codeletsPerTP2409), TP2558Ptr(new TP2558 *[NUMTPS2558]), TP2558_alreadyLaunched(new size_t [NUMTPS2558]), numTPsSet2558(0), numTPsReady2558(0), TPsToUse2558(NUMTPS2558), codeletsPerTP2558(this->numThreads/NUMTPS2558), totalCodelets2558(this->TPsToUse2558*this->codeletsPerTP2558), TP3196Ptr(new TP3196 *[NUMTPS3196]), TP3196_alreadyLaunched(new size_t [NUMTPS3196]), numTPsSet3196(0), numTPsReady3196(0), TPsToUse3196(NUMTPS3196), codeletsPerTP3196(this->numThreads/NUMTPS3196), totalCodelets3196(this->TPsToUse3196*this->codeletsPerTP3196), TP3345Ptr(new TP3345 *[NUMTPS3345]), TP3345_alreadyLaunched(new size_t [NUMTPS3345]), numTPsSet3345(0), numTPsReady3345(0), TPsToUse3345(NUMTPS3345), codeletsPerTP3345(this->numThreads/NUMTPS3345), totalCodelets3345(this->TPsToUse3345*this->codeletsPerTP3345), TP3980Ptr(new TP3980 *[NUMTPS3980]), TP3980_alreadyLaunched(new size_t [NUMTPS3980]), numTPsSet3980(0), numTPsReady3980(0), TPsToUse3980(NUMTPS3980), codeletsPerTP3980(this->numThreads/NUMTPS3980), totalCodelets3980(this->TPsToUse3980*this->codeletsPerTP3980) ,barrierCodelets2204(new _barrierCodelets2204[1]) ,checkInCodelets2206(new _checkInCodelets2206[this->numThreads]) ,checkInCodelets2223(new _checkInCodelets2223[this->numThreads]) ,barrierCodelets2223(new _barrierCodelets2223[1]) ,checkInCodelets2274(new _checkInCodelets2274[this->numThreads]) ,barrierCodelets2274(new _barrierCodelets2274[1]) ,checkInCodelets2406(new _checkInCodelets2406[this->numThreads]) ,checkInCodelets2409(new _checkInCodelets2409[this->numThreads]) ,barrierCodelets2409(new _barrierCodelets2409[1]) ,checkInCodelets2558(new _checkInCodelets2558[this->numThreads]) ,barrierCodelets2558(new _barrierCodelets2558[1]) ,checkInCodelets3193(new _checkInCodelets3193[this->numThreads]) ,checkInCodelets3196(new _checkInCodelets3196[this->numThreads]) ,barrierCodelets3196(new _barrierCodelets3196[1]) ,checkInCodelets3345(new _checkInCodelets3345[this->numThreads]) ,barrierCodelets3345(new _barrierCodelets3345[1]) ,checkInCodelets3980(new _checkInCodelets3980[this->numThreads]) ,barrierCodelets3980(new _barrierCodelets3980[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets2204[0] = _barrierCodelets2204(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets3980[0] = _barrierCodelets3980(NUMTPS3980,NUMTPS3980,this, 0);
barrierCodelets3345[0] = _barrierCodelets3345(NUMTPS3345,NUMTPS3345,this, 0);
barrierCodelets3196[0] = _barrierCodelets3196(NUMTPS3196,NUMTPS3196,this, 0);
barrierCodelets2558[0] = _barrierCodelets2558(NUMTPS2558,NUMTPS2558,this, 0);
barrierCodelets2409[0] = _barrierCodelets2409(NUMTPS2409,NUMTPS2409,this, 0);
barrierCodelets2274[0] = _barrierCodelets2274(NUMTPS2274,NUMTPS2274,this, 0);
barrierCodelets2223[0] = _barrierCodelets2223(NUMTPS2223,NUMTPS2223,this, 0);
_checkInCodelets3980 * checkInCodelets3980Ptr = (this->checkInCodelets3980);
for(int i=0; i<NUMTPS3980; i++)
{
TP3980Ptr[i] = nullptr;
TP3980_alreadyLaunched[i] = 0;
}
_checkInCodelets3345 * checkInCodelets3345Ptr = (this->checkInCodelets3345);
for(int i=0; i<NUMTPS3345; i++)
{
TP3345Ptr[i] = nullptr;
TP3345_alreadyLaunched[i] = 0;
}
_checkInCodelets3196 * checkInCodelets3196Ptr = (this->checkInCodelets3196);
for(int i=0; i<NUMTPS3196; i++)
{
TP3196Ptr[i] = nullptr;
TP3196_alreadyLaunched[i] = 0;
}
_checkInCodelets3193 * checkInCodelets3193Ptr = (this->checkInCodelets3193);
_checkInCodelets2558 * checkInCodelets2558Ptr = (this->checkInCodelets2558);
for(int i=0; i<NUMTPS2558; i++)
{
TP2558Ptr[i] = nullptr;
TP2558_alreadyLaunched[i] = 0;
}
_checkInCodelets2409 * checkInCodelets2409Ptr = (this->checkInCodelets2409);
for(int i=0; i<NUMTPS2409; i++)
{
TP2409Ptr[i] = nullptr;
TP2409_alreadyLaunched[i] = 0;
}
_checkInCodelets2406 * checkInCodelets2406Ptr = (this->checkInCodelets2406);
_checkInCodelets2274 * checkInCodelets2274Ptr = (this->checkInCodelets2274);
for(int i=0; i<NUMTPS2274; i++)
{
TP2274Ptr[i] = nullptr;
TP2274_alreadyLaunched[i] = 0;
}
_checkInCodelets2223 * checkInCodelets2223Ptr = (this->checkInCodelets2223);
for(int i=0; i<NUMTPS2223; i++)
{
TP2223Ptr[i] = nullptr;
TP2223_alreadyLaunched[i] = 0;
}
_checkInCodelets2206 * checkInCodelets2206Ptr = (this->checkInCodelets2206);
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets3980Ptr) = _checkInCodelets3980(1,1,this,codeletCounter);
checkInCodelets3980Ptr++;
(*checkInCodelets3345Ptr) = _checkInCodelets3345(1,1,this,codeletCounter);
checkInCodelets3345Ptr++;
(*checkInCodelets3196Ptr) = _checkInCodelets3196(1,1,this,codeletCounter);
checkInCodelets3196Ptr++;
(*checkInCodelets3193Ptr) = _checkInCodelets3193(1,1,this,codeletCounter);
checkInCodelets3193Ptr++;
(*checkInCodelets2558Ptr) = _checkInCodelets2558(1,1,this,codeletCounter);
checkInCodelets2558Ptr++;
(*checkInCodelets2409Ptr) = _checkInCodelets2409(1,1,this,codeletCounter);
checkInCodelets2409Ptr++;
(*checkInCodelets2406Ptr) = _checkInCodelets2406(1,1,this,codeletCounter);
checkInCodelets2406Ptr++;
(*checkInCodelets2274Ptr) = _checkInCodelets2274(1,1,this,codeletCounter);
checkInCodelets2274Ptr++;
(*checkInCodelets2223Ptr) = _checkInCodelets2223(1,1,this,codeletCounter);
checkInCodelets2223Ptr++;
(*checkInCodelets2206Ptr) = _checkInCodelets2206(1,1,this,codeletCounter);
(*checkInCodelets2206Ptr).decDep();
checkInCodelets2206Ptr++;
}
}
TP2204::~TP2204(){
delete []L1_darts2204;
delete []L2_darts2204;
delete []dsspm_darts2204;
delete []eta_darts2204;
delete []i_darts2204;
delete []iend1_darts2204;
delete []iglob_darts2204;
delete []ist1_darts2204;
delete []j_darts2204;
delete []jend1_darts2204;
delete []jglob_darts2204;
delete []jst1_darts2204;
delete []k_darts2204;
delete []m_darts2204;
delete []q_darts2204;
delete []tmp_darts2204;
delete []u21_darts2204;
delete []u21i_darts2204;
delete []u21im1_darts2204;
delete []u21j_darts2204;
delete []u21jm1_darts2204;
delete []u21k_darts2204;
delete []u21km1_darts2204;
delete []u31_darts2204;
delete []u31i_darts2204;
delete []u31im1_darts2204;
delete []u31j_darts2204;
delete []u31jm1_darts2204;
delete []u31k_darts2204;
delete []u31km1_darts2204;
delete []u41_darts2204;
delete []u41i_darts2204;
delete []u41im1_darts2204;
delete []u41j_darts2204;
delete []u41jm1_darts2204;
delete []u41k_darts2204;
delete []u41km1_darts2204;
delete []u51i_darts2204;
delete []u51im1_darts2204;
delete []u51j_darts2204;
delete []u51jm1_darts2204;
delete []u51k_darts2204;
delete []u51km1_darts2204;
delete []xi_darts2204;
delete []zeta_darts2204;
delete [] barrierCodelets2204;
delete [] barrierCodelets3980;
delete [] checkInCodelets3980;
delete [] barrierCodelets3345;
delete [] checkInCodelets3345;
delete [] barrierCodelets3196;
delete [] checkInCodelets3196;
delete [] checkInCodelets3193;
delete [] barrierCodelets2558;
delete [] checkInCodelets2558;
delete [] barrierCodelets2409;
delete [] checkInCodelets2409;
delete [] checkInCodelets2406;
delete [] barrierCodelets2274;
delete [] checkInCodelets2274;
delete [] barrierCodelets2223;
delete [] checkInCodelets2223;
delete [] checkInCodelets2206;
}
/*TP2223: OMPForDirective*/
void TP2223::_barrierCodelets2223::fire(void)
{
TP2223* myTP = static_cast<TP2223*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets2223[0].decDep ();
}
bool TP2223::requestNewRangeIterations2223(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet2223 * codeletID;
int tempEndRange   = rangePerCodelet2223 * (codeletID + 1);
if (remainderRange2223 != 0)
{
if (codeletID < (uint32_t)remainderRange2223)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange2223;
tempEndRange += remainderRange2223;
}
}
tempStartRange = tempStartRange*1 + minIteration2223;
tempEndRange = tempEndRange*1 + minIteration2223;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration2223 < lastIteration2223)
{
(this->inputsTPParent->i_darts2223[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts2223[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration2223;
}
}
return isThereNewIteration;
}
void TP2223::_checkInCodelets2224::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 2224: ForStmt*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: m*/
int* i = &(this->inputsTPParent->i_darts2223[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts2223[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts2223[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts2223[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2223((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2223[0].decDep();
return;
}
for (int i_darts_counter_temp2223 = (*i);i_darts_counter_temp2223<endRange && i_darts_counter_temp2223<this->inputsTPParent->lastIteration2223;i_darts_counter_temp2223++)
{
{
{
/*Loop's init*/
(*j) = 0;
int j_darts_counter_temp2223 = (*j);
for(;j_darts_counter_temp2223 < ny;j_darts_counter_temp2223++){
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp2223 = (*k);
for(;k_darts_counter_temp2223 < nz;k_darts_counter_temp2223++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2223 = (*m);
for(;m_darts_counter_temp2223 < 5;m_darts_counter_temp2223++){
frct[(i_darts_counter_temp2223)][j_darts_counter_temp2223][k_darts_counter_temp2223][m_darts_counter_temp2223] = 0.;
}
(*m) = m_darts_counter_temp2223;
}
}
(*k) = k_darts_counter_temp2223;
}
}
(*j) = j_darts_counter_temp2223;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2223[0].decDep();
}
TP2223::TP2223(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2223** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts2223(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts2223(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts2223(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts2223(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration2223(in_initIteration), lastIteration2223(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets2223(new _barrierCodelets2223[1]) ,checkInCodelets2224(new _checkInCodelets2224[this->numThreads]){
/*Initialize the loop parameters*/
range2223 = abs (lastIteration2223 - initIteration2223) / 1;
rangePerCodelet2223 = range2223 / numThreads;
minIteration2223 = min<int>(lastIteration2223, initIteration2223);
remainderRange2223 = range2223 % numThreads;
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets2223[0] = _barrierCodelets2223(this->numThreads,this->numThreads,this, 0);
_checkInCodelets2224 * checkInCodelets2224Ptr = (this->checkInCodelets2224);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets2224);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets2224Ptr) = _checkInCodelets2224(2,1,this,codeletCounter);
#else
(*checkInCodelets2224Ptr) = _checkInCodelets2224(1,1,this,codeletCounter);
#endif
(*checkInCodelets2224Ptr).decDep();
checkInCodelets2224Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP2223::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets2224[localID].setID (codeletID);
this->checkInCodelets2224[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets2224[localID + this->baseNumThreads * i] = _checkInCodelets2224(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets2224[localID + this->baseNumThreads * i] = _checkInCodelets2224(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets2224[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets2224[localID + this->baseNumThreads * i].decDep();
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
TP2223::~TP2223(){
delete [] barrierCodelets2223;
delete [] checkInCodelets2224;
}
/*TP2274: OMPForDirective*/
void TP2274::_barrierCodelets2274::fire(void)
{
TP2274* myTP = static_cast<TP2274*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets2274[0].decDep ();
}
bool TP2274::requestNewRangeIterations2274(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet2274 * codeletID;
int tempEndRange   = rangePerCodelet2274 * (codeletID + 1);
if (remainderRange2274 != 0)
{
if (codeletID < (uint32_t)remainderRange2274)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange2274;
tempEndRange += remainderRange2274;
}
}
tempStartRange = tempStartRange*1 + minIteration2274;
tempEndRange = tempEndRange*1 + minIteration2274;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration2274 < lastIteration2274)
{
(this->inputsTPParent->i_darts2274[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts2274[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration2274;
}
}
return isThereNewIteration;
}
void TP2274::_checkInCodelets2275::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->eta_darts2274[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->eta_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iglob_darts2274[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jglob_darts2274[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->xi_darts2274[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->xi_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->zeta_darts2274[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts2204[this->getID()]);

/*printing node 2275: ForStmt*/
/*var: eta*/
/*var: i*/
/*var: iglob*/
/*var: j*/
/*var: jglob*/
/*var: k*/
/*var: m*/
/*var: xi*/
/*var: zeta*/
double** eta = &(this->inputsTPParent->eta_darts2274[this->getLocalID()]);
(void)eta/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts2274[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iglob = &(this->inputsTPParent->iglob_darts2274[this->getLocalID()]);
(void)iglob/*OMP_SHARED_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts2274[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jglob = &(this->inputsTPParent->jglob_darts2274[this->getLocalID()]);
(void)jglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts2274[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts2274[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** xi = &(this->inputsTPParent->xi_darts2274[this->getLocalID()]);
(void)xi/*OMP_SHARED_PRIVATE*/;
double** zeta = &(this->inputsTPParent->zeta_darts2274[this->getLocalID()]);
(void)zeta/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2274((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2274[0].decDep();
return;
}
for (int i_darts_counter_temp2274 = (*i);i_darts_counter_temp2274<endRange && i_darts_counter_temp2274<this->inputsTPParent->lastIteration2274;i_darts_counter_temp2274++)
{
{
(*(*iglob)) = (i_darts_counter_temp2274);
(*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
{
/*Loop's init*/
(*j) = 0;
int j_darts_counter_temp2274 = (*j);
for(;j_darts_counter_temp2274 < ny;j_darts_counter_temp2274++){
(*(*jglob)) = j_darts_counter_temp2274;
(*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp2274 = (*k);
for(;k_darts_counter_temp2274 < nz;k_darts_counter_temp2274++){
(*(*zeta)) = ((double)(k_darts_counter_temp2274)) / (nz - 1);
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2274 = (*m);
for(;m_darts_counter_temp2274 < 5;m_darts_counter_temp2274++){
rsd[(i_darts_counter_temp2274)][j_darts_counter_temp2274][k_darts_counter_temp2274][m_darts_counter_temp2274] = ce[m_darts_counter_temp2274][0] + ce[m_darts_counter_temp2274][1] * (*(*xi)) + ce[m_darts_counter_temp2274][2] * (*(*eta)) + ce[m_darts_counter_temp2274][3] * (*(*zeta)) + ce[m_darts_counter_temp2274][4] * (*(*xi)) * (*(*xi)) + ce[m_darts_counter_temp2274][5] * (*(*eta)) * (*(*eta)) + ce[m_darts_counter_temp2274][6] * (*(*zeta)) * (*(*zeta)) + ce[m_darts_counter_temp2274][7] * (*(*xi)) * (*(*xi)) * (*(*xi)) + ce[m_darts_counter_temp2274][8] * (*(*eta)) * (*(*eta)) * (*(*eta)) + ce[m_darts_counter_temp2274][9] * (*(*zeta)) * (*(*zeta)) * (*(*zeta)) + ce[m_darts_counter_temp2274][10] * (*(*xi)) * (*(*xi)) * (*(*xi)) * (*(*xi)) + ce[m_darts_counter_temp2274][11] * (*(*eta)) * (*(*eta)) * (*(*eta)) * (*(*eta)) + ce[m_darts_counter_temp2274][12] * (*(*zeta)) * (*(*zeta)) * (*(*zeta)) * (*(*zeta));
}
(*m) = m_darts_counter_temp2274;
}
}
(*k) = k_darts_counter_temp2274;
}
}
(*j) = j_darts_counter_temp2274;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2274[0].decDep();
}
TP2274::TP2274(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2274** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),eta_darts2274(new double*[this->numThreads]),i_darts2274(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iglob_darts2274(new int*[this->numThreads]),j_darts2274(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jglob_darts2274(new int*[this->numThreads]),k_darts2274(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts2274(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,xi_darts2274(new double*[this->numThreads]),zeta_darts2274(new double*[this->numThreads]), initIteration2274(in_initIteration), lastIteration2274(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets2274(new _barrierCodelets2274[1]) ,checkInCodelets2275(new _checkInCodelets2275[this->numThreads]){
/*Initialize the loop parameters*/
range2274 = abs (lastIteration2274 - initIteration2274) / 1;
rangePerCodelet2274 = range2274 / numThreads;
minIteration2274 = min<int>(lastIteration2274, initIteration2274);
remainderRange2274 = range2274 % numThreads;
/*Initialize inputs and vars.*/
this->eta_darts2274 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->iglob_darts2274 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jglob_darts2274 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->xi_darts2274 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->zeta_darts2274 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets2274[0] = _barrierCodelets2274(this->numThreads,this->numThreads,this, 0);
_checkInCodelets2275 * checkInCodelets2275Ptr = (this->checkInCodelets2275);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets2275);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets2275Ptr) = _checkInCodelets2275(2,1,this,codeletCounter);
#else
(*checkInCodelets2275Ptr) = _checkInCodelets2275(1,1,this,codeletCounter);
#endif
(*checkInCodelets2275Ptr).decDep();
checkInCodelets2275Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP2274::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets2275[localID].setID (codeletID);
this->checkInCodelets2275[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets2275[localID + this->baseNumThreads * i] = _checkInCodelets2275(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets2275[localID + this->baseNumThreads * i] = _checkInCodelets2275(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets2275[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets2275[localID + this->baseNumThreads * i].decDep();
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
TP2274::~TP2274(){
delete [] eta_darts2274;
delete [] iglob_darts2274;
delete [] jglob_darts2274;
delete [] xi_darts2274;
delete [] zeta_darts2274;
delete [] barrierCodelets2274;
delete [] checkInCodelets2275;
}
/*TP2409: OMPForDirective*/
void TP2409::_barrierCodelets2409::fire(void)
{
TP2409* myTP = static_cast<TP2409*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets2409[0].decDep ();
}
bool TP2409::requestNewRangeIterations2409(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet2409 * codeletID;
int tempEndRange   = rangePerCodelet2409 * (codeletID + 1);
if (remainderRange2409 != 0)
{
if (codeletID < (uint32_t)remainderRange2409)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange2409;
tempEndRange += remainderRange2409;
}
}
tempStartRange = tempStartRange*1 + minIteration2409;
tempEndRange = tempEndRange*1 + minIteration2409;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration2409 < lastIteration2409)
{
(this->inputsTPParent->i_darts2409[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts2409[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration2409;
}
}
return isThereNewIteration;
}
void TP2409::_checkInCodelets2410::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L1_darts2409[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts2409[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts2409[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21_darts2409[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21_darts2204[this->getID()]);

/*printing node 2410: ForStmt*/
/*var: L1*/
/*var: L2*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: q*/
/*var: u21*/
int* i = &(this->inputsTPParent->i_darts2409[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts2409[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts2409[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts2409[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** u21 = &(this->inputsTPParent->u21_darts2409[this->getLocalID()]);
(void)u21/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2409((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2409[0].decDep();
return;
}
for (int i_darts_counter_temp2409 = (*i);i_darts_counter_temp2409<=endRange && i_darts_counter_temp2409<=this->inputsTPParent->lastIteration2409;i_darts_counter_temp2409++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp2409 = (*j);
for(;j_darts_counter_temp2409 <= jend;j_darts_counter_temp2409++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp2409 = (*k);
for(;k_darts_counter_temp2409 < nz - 1;k_darts_counter_temp2409++){
flux[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][0] = rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1];
(*(*u21)) = rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1] / rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][0];
(*(*q)) = 0.5 * (rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1] * rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1] + rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][2] * rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][2] + rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][3] * rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][3]) / rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][0];
flux[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1] = rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][1] * (*(*u21)) + 0.40000000000000002 * (rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][4] - (*(*q)));
flux[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][2] = rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][2] * (*(*u21));
flux[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][3] = rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][3] * (*(*u21));
flux[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][4] = (1.3999999999999999 * rsd[(i_darts_counter_temp2409)][j_darts_counter_temp2409][k_darts_counter_temp2409][4] - 0.40000000000000002 * (*(*q))) * (*(*u21));
}
(*k) = k_darts_counter_temp2409;
}
}
(*j) = j_darts_counter_temp2409;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2409[0].decDep();
}
TP2409::TP2409(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2409** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L1_darts2409(new int*[this->numThreads]),L2_darts2409(new int*[this->numThreads]),i_darts2409(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts2409(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts2409(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts2409(new double*[this->numThreads]),u21_darts2409(new double*[this->numThreads]), initIteration2409(in_initIteration), lastIteration2409(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets2409(new _barrierCodelets2409[1]) ,checkInCodelets2410(new _checkInCodelets2410[this->numThreads]){
/*Initialize the loop parameters*/
range2409 = abs (lastIteration2409 - initIteration2409) / 1;
rangePerCodelet2409 = range2409 / numThreads;
minIteration2409 = min<int>(lastIteration2409, initIteration2409);
remainderRange2409 = range2409 % numThreads;
/*Initialize inputs and vars.*/
this->L1_darts2409 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->L2_darts2409 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->q_darts2409 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21_darts2409 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets2409[0] = _barrierCodelets2409(this->numThreads,this->numThreads,this, 0);
_checkInCodelets2410 * checkInCodelets2410Ptr = (this->checkInCodelets2410);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets2410);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets2410Ptr) = _checkInCodelets2410(2,1,this,codeletCounter);
#else
(*checkInCodelets2410Ptr) = _checkInCodelets2410(1,1,this,codeletCounter);
#endif
(*checkInCodelets2410Ptr).decDep();
checkInCodelets2410Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP2409::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets2410[localID].setID (codeletID);
this->checkInCodelets2410[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets2410[localID + this->baseNumThreads * i] = _checkInCodelets2410(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets2410[localID + this->baseNumThreads * i] = _checkInCodelets2410(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets2410[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets2410[localID + this->baseNumThreads * i].decDep();
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
TP2409::~TP2409(){
delete [] L1_darts2409;
delete [] L2_darts2409;
delete [] q_darts2409;
delete [] u21_darts2409;
delete [] barrierCodelets2409;
delete [] checkInCodelets2410;
}
/*TP2558: OMPForDirective*/
void TP2558::_barrierCodelets2558::fire(void)
{
TP2558* myTP = static_cast<TP2558*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets2558[0].decDep ();
}
bool TP2558::requestNewRangeIterations2558(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet2558 * codeletID;
int tempEndRange   = rangePerCodelet2558 * (codeletID + 1);
if (remainderRange2558 != 0)
{
if (codeletID < (uint32_t)remainderRange2558)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange2558;
tempEndRange += remainderRange2558;
}
}
tempStartRange = tempStartRange*1 + minIteration2558;
tempEndRange = tempEndRange*1 + minIteration2558;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration2558 < lastIteration2558)
{
(this->inputsTPParent->j_darts2558[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->j_darts2558[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration2558;
}
}
return isThereNewIteration;
}
void TP2558::_checkInCodelets2559::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts2558[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->dsspm_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iend1_darts2558[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->ist1_darts2558[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21i_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21im1_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31i_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31im1_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41i_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41im1_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51i_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51im1_darts2558[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts2204[this->getID()]);

/*printing node 2559: ForStmt*/
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
int** L2 = &(this->inputsTPParent->L2_darts2558[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
double** dsspm = &(this->inputsTPParent->dsspm_darts2558[this->getLocalID()]);
(void)dsspm/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts2558[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iend1 = &(this->inputsTPParent->iend1_darts2558[this->getLocalID()]);
(void)iend1/*OMP_SHARED_PRIVATE*/;
int** ist1 = &(this->inputsTPParent->ist1_darts2558[this->getLocalID()]);
(void)ist1/*OMP_SHARED_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts2558[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts2558[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts2558[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts2558[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21i = &(this->inputsTPParent->u21i_darts2558[this->getLocalID()]);
(void)u21i/*OMP_SHARED_PRIVATE*/;
double** u21im1 = &(this->inputsTPParent->u21im1_darts2558[this->getLocalID()]);
(void)u21im1/*OMP_SHARED_PRIVATE*/;
double** u31i = &(this->inputsTPParent->u31i_darts2558[this->getLocalID()]);
(void)u31i/*OMP_SHARED_PRIVATE*/;
double** u31im1 = &(this->inputsTPParent->u31im1_darts2558[this->getLocalID()]);
(void)u31im1/*OMP_SHARED_PRIVATE*/;
double** u41i = &(this->inputsTPParent->u41i_darts2558[this->getLocalID()]);
(void)u41i/*OMP_SHARED_PRIVATE*/;
double** u41im1 = &(this->inputsTPParent->u41im1_darts2558[this->getLocalID()]);
(void)u41im1/*OMP_SHARED_PRIVATE*/;
double** u51i = &(this->inputsTPParent->u51i_darts2558[this->getLocalID()]);
(void)u51i/*OMP_SHARED_PRIVATE*/;
double** u51im1 = &(this->inputsTPParent->u51im1_darts2558[this->getLocalID()]);
(void)u51im1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2558((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2558[0].decDep();
return;
}
for (int j_darts_counter_temp2558 = (*j);j_darts_counter_temp2558<=endRange && j_darts_counter_temp2558<=this->inputsTPParent->lastIteration2558;j_darts_counter_temp2558++)
{
{
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp2558 = (*k);
for(;k_darts_counter_temp2558 <= nz - 2;k_darts_counter_temp2558++){
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp2558 = (*i);
for(;i_darts_counter_temp2558 <= iend;i_darts_counter_temp2558++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2558 = (*m);
for(;m_darts_counter_temp2558 < 5;m_darts_counter_temp2558++){
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - tx2 * (flux[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - flux[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
}
(*m) = m_darts_counter_temp2558;
}
}
(*i) = i_darts_counter_temp2558;
}
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp2558 = (*i);
for(;i_darts_counter_temp2558 <= (*(*L2));i_darts_counter_temp2558++){
(*(*tmp)) = 1. / rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0];
(*(*u21i)) = (*(*tmp)) * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1];
(*(*u31i)) = (*(*tmp)) * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2];
(*(*u41i)) = (*(*tmp)) * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3];
(*(*u51i)) = (*(*tmp)) * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4];
(*(*tmp)) = 1. / rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0];
(*(*u21im1)) = (*(*tmp)) * rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1];
(*(*u31im1)) = (*(*tmp)) * rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2];
(*(*u41im1)) = (*(*tmp)) * rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3];
(*(*u51im1)) = (*(*tmp)) * rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4];
flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] = tx3 * ((*(*u31i)) - (*(*u31im1)));
flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] = tx3 * ((*(*u41i)) - (*(*u41im1)));
flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tx3 * (((*(*u21i)) * (*(*u21i)) + (*(*u31i)) * (*(*u31i)) + (*(*u41i)) * (*(*u41i))) - ((*(*u21im1)) * (*(*u21im1)) + (*(*u31im1)) * (*(*u31im1)) + (*(*u41im1)) * (*(*u41im1)))) + (1. / 6.) * tx3 * ((*(*u21i)) * (*(*u21i)) - (*(*u21im1)) * (*(*u21im1))) + 1.3999999999999999 * 1.3999999999999999 * tx3 * ((*(*u51i)) - (*(*u51im1)));
}
(*i) = i_darts_counter_temp2558;
}
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp2558 = (*i);
for(;i_darts_counter_temp2558 <= iend;i_darts_counter_temp2558++){
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0] + dx1 * tx1 * (rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0] - 2. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0] + rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][0]);
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] - flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1]) + dx2 * tx1 * (rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] - 2. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1] + rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][1]);
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] - flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2]) + dx3 * tx1 * (rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] - 2. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2] + rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][2]);
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] - flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3]) + dx4 * tx1 * (rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] - 2. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3] + rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][3]);
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] - flux[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4]) + dx5 * tx1 * (rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] - 2. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4] + rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][4]);
}
(*i) = i_darts_counter_temp2558;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2558 = (*m);
for(;m_darts_counter_temp2558 < 5;m_darts_counter_temp2558++){
frct[1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - (*(*dsspm)) * (+5. * rsd[1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + rsd[3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
frct[2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - (*(*dsspm)) * (-4. * rsd[1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + 6. * rsd[2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + rsd[4][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
}
(*m) = m_darts_counter_temp2558;
}
(*(*ist1)) = 3;
(*(*iend1)) = nx - 4;
{
/*Loop's init*/
(*i) = (*(*ist1));
int i_darts_counter_temp2558 = (*i);
for(;i_darts_counter_temp2558 <= (*(*iend1));i_darts_counter_temp2558++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2558 = (*m);
for(;m_darts_counter_temp2558 < 5;m_darts_counter_temp2558++){
frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - (*(*dsspm)) * (rsd[i_darts_counter_temp2558 - 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[i_darts_counter_temp2558 - 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + 6. * rsd[i_darts_counter_temp2558][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[i_darts_counter_temp2558 + 1][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + rsd[i_darts_counter_temp2558 + 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
}
(*m) = m_darts_counter_temp2558;
}
}
(*i) = i_darts_counter_temp2558;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp2558 = (*m);
for(;m_darts_counter_temp2558 < 5;m_darts_counter_temp2558++){
frct[nx - 3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[nx - 3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - (*(*dsspm)) * (rsd[nx - 5][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[nx - 4][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + 6. * rsd[nx - 3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[nx - 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
frct[nx - 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] = frct[nx - 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - (*(*dsspm)) * (rsd[nx - 4][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] - 4. * rsd[nx - 3][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558] + 5. * rsd[nx - 2][(j_darts_counter_temp2558)][k_darts_counter_temp2558][m_darts_counter_temp2558]);
}
(*m) = m_darts_counter_temp2558;
}
}
(*k) = k_darts_counter_temp2558;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets2558[0].decDep();
}
TP2558::TP2558(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2558** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L2_darts2558(new int*[this->numThreads]),dsspm_darts2558(new double*[this->numThreads]),i_darts2558(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iend1_darts2558(new int*[this->numThreads]),ist1_darts2558(new int*[this->numThreads]),j_darts2558(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts2558(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts2558(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts2558(new double*[this->numThreads]),u21i_darts2558(new double*[this->numThreads]),u21im1_darts2558(new double*[this->numThreads]),u31i_darts2558(new double*[this->numThreads]),u31im1_darts2558(new double*[this->numThreads]),u41i_darts2558(new double*[this->numThreads]),u41im1_darts2558(new double*[this->numThreads]),u51i_darts2558(new double*[this->numThreads]),u51im1_darts2558(new double*[this->numThreads]), initIteration2558(in_initIteration), lastIteration2558(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets2558(new _barrierCodelets2558[1]) ,checkInCodelets2559(new _checkInCodelets2559[this->numThreads]){
/*Initialize the loop parameters*/
range2558 = abs (lastIteration2558 - initIteration2558) / 1;
rangePerCodelet2558 = range2558 / numThreads;
minIteration2558 = min<int>(lastIteration2558, initIteration2558);
remainderRange2558 = range2558 % numThreads;
/*Initialize inputs and vars.*/
this->L2_darts2558 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->dsspm_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->iend1_darts2558 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->ist1_darts2558 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21i_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21im1_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31i_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31im1_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41i_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41im1_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51i_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51im1_darts2558 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets2558[0] = _barrierCodelets2558(this->numThreads,this->numThreads,this, 0);
_checkInCodelets2559 * checkInCodelets2559Ptr = (this->checkInCodelets2559);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets2559);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets2559Ptr) = _checkInCodelets2559(2,1,this,codeletCounter);
#else
(*checkInCodelets2559Ptr) = _checkInCodelets2559(1,1,this,codeletCounter);
#endif
(*checkInCodelets2559Ptr).decDep();
checkInCodelets2559Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP2558::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets2559[localID].setID (codeletID);
this->checkInCodelets2559[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets2559[localID + this->baseNumThreads * i] = _checkInCodelets2559(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets2559[localID + this->baseNumThreads * i] = _checkInCodelets2559(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets2559[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets2559[localID + this->baseNumThreads * i].decDep();
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
TP2558::~TP2558(){
delete [] L2_darts2558;
delete [] dsspm_darts2558;
delete [] iend1_darts2558;
delete [] ist1_darts2558;
delete [] tmp_darts2558;
delete [] u21i_darts2558;
delete [] u21im1_darts2558;
delete [] u31i_darts2558;
delete [] u31im1_darts2558;
delete [] u41i_darts2558;
delete [] u41im1_darts2558;
delete [] u51i_darts2558;
delete [] u51im1_darts2558;
delete [] barrierCodelets2558;
delete [] checkInCodelets2559;
}
/*TP3196: OMPForDirective*/
void TP3196::_barrierCodelets3196::fire(void)
{
TP3196* myTP = static_cast<TP3196*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets3196[0].decDep ();
}
bool TP3196::requestNewRangeIterations3196(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet3196 * codeletID;
int tempEndRange   = rangePerCodelet3196 * (codeletID + 1);
if (remainderRange3196 != 0)
{
if (codeletID < (uint32_t)remainderRange3196)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange3196;
tempEndRange += remainderRange3196;
}
}
tempStartRange = tempStartRange*1 + minIteration3196;
tempEndRange = tempEndRange*1 + minIteration3196;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration3196 < lastIteration3196)
{
(this->inputsTPParent->i_darts3196[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts3196[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration3196;
}
}
return isThereNewIteration;
}
void TP3196::_checkInCodelets3197::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L1_darts3196[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts3196[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts3196[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31_darts3196[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31_darts2204[this->getID()]);

/*printing node 3197: ForStmt*/
/*var: L1*/
/*var: L2*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: q*/
/*var: u31*/
int** L1 = &(this->inputsTPParent->L1_darts3196[this->getLocalID()]);
(void)L1/*OMP_SHARED_PRIVATE*/;
int** L2 = &(this->inputsTPParent->L2_darts3196[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts3196[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts3196[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts3196[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts3196[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** u31 = &(this->inputsTPParent->u31_darts3196[this->getLocalID()]);
(void)u31/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3196((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3196[0].decDep();
return;
}
for (int i_darts_counter_temp3196 = (*i);i_darts_counter_temp3196<=endRange && i_darts_counter_temp3196<=this->inputsTPParent->lastIteration3196;i_darts_counter_temp3196++)
{
{
{
/*Loop's init*/
(*j) = (*(*L1));
int j_darts_counter_temp3196 = (*j);
for(;j_darts_counter_temp3196 <= (*(*L2));j_darts_counter_temp3196++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp3196 = (*k);
for(;k_darts_counter_temp3196 <= nz - 2;k_darts_counter_temp3196++){
flux[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][0] = rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2];
(*(*u31)) = rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2] / rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][0];
(*(*q)) = 0.5 * (rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][1] * rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][1] + rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2] * rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2] + rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][3] * rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][3]) / rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][0];
flux[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][1] = rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][1] * (*(*u31));
flux[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2] = rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][2] * (*(*u31)) + 0.40000000000000002 * (rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][4] - (*(*q)));
flux[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][3] = rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][3] * (*(*u31));
flux[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][4] = (1.3999999999999999 * rsd[(i_darts_counter_temp3196)][j_darts_counter_temp3196][k_darts_counter_temp3196][4] - 0.40000000000000002 * (*(*q))) * (*(*u31));
}
(*k) = k_darts_counter_temp3196;
}
}
(*j) = j_darts_counter_temp3196;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3196[0].decDep();
}
TP3196::TP3196(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3196** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L1_darts3196(new int*[this->numThreads]),L2_darts3196(new int*[this->numThreads]),i_darts3196(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts3196(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts3196(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts3196(new double*[this->numThreads]),u31_darts3196(new double*[this->numThreads]), initIteration3196(in_initIteration), lastIteration3196(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets3196(new _barrierCodelets3196[1]) ,checkInCodelets3197(new _checkInCodelets3197[this->numThreads]){
/*Initialize the loop parameters*/
range3196 = abs (lastIteration3196 - initIteration3196) / 1;
rangePerCodelet3196 = range3196 / numThreads;
minIteration3196 = min<int>(lastIteration3196, initIteration3196);
remainderRange3196 = range3196 % numThreads;
/*Initialize inputs and vars.*/
this->L1_darts3196 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->L2_darts3196 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->q_darts3196 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31_darts3196 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets3196[0] = _barrierCodelets3196(this->numThreads,this->numThreads,this, 0);
_checkInCodelets3197 * checkInCodelets3197Ptr = (this->checkInCodelets3197);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets3197);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets3197Ptr) = _checkInCodelets3197(2,1,this,codeletCounter);
#else
(*checkInCodelets3197Ptr) = _checkInCodelets3197(1,1,this,codeletCounter);
#endif
(*checkInCodelets3197Ptr).decDep();
checkInCodelets3197Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP3196::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets3197[localID].setID (codeletID);
this->checkInCodelets3197[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets3197[localID + this->baseNumThreads * i] = _checkInCodelets3197(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets3197[localID + this->baseNumThreads * i] = _checkInCodelets3197(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets3197[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets3197[localID + this->baseNumThreads * i].decDep();
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
TP3196::~TP3196(){
delete [] L1_darts3196;
delete [] L2_darts3196;
delete [] q_darts3196;
delete [] u31_darts3196;
delete [] barrierCodelets3196;
delete [] checkInCodelets3197;
}
/*TP3345: OMPForDirective*/
void TP3345::_barrierCodelets3345::fire(void)
{
TP3345* myTP = static_cast<TP3345*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets3345[0].decDep ();
}
bool TP3345::requestNewRangeIterations3345(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet3345 * codeletID;
int tempEndRange   = rangePerCodelet3345 * (codeletID + 1);
if (remainderRange3345 != 0)
{
if (codeletID < (uint32_t)remainderRange3345)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange3345;
tempEndRange += remainderRange3345;
}
}
tempStartRange = tempStartRange*1 + minIteration3345;
tempEndRange = tempEndRange*1 + minIteration3345;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration3345 < lastIteration3345)
{
(this->inputsTPParent->i_darts3345[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts3345[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration3345;
}
}
return isThereNewIteration;
}
void TP3345::_checkInCodelets3346::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts3345[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->dsspm_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jend1_darts3345[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jst1_darts3345[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21j_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21jm1_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31j_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31jm1_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41j_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41jm1_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51j_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51jm1_darts3345[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts2204[this->getID()]);

/*printing node 3346: ForStmt*/
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
int** L2 = &(this->inputsTPParent->L2_darts3345[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
double** dsspm = &(this->inputsTPParent->dsspm_darts3345[this->getLocalID()]);
(void)dsspm/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts3345[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts3345[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jend1 = &(this->inputsTPParent->jend1_darts3345[this->getLocalID()]);
(void)jend1/*OMP_SHARED_PRIVATE*/;
int** jst1 = &(this->inputsTPParent->jst1_darts3345[this->getLocalID()]);
(void)jst1/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts3345[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts3345[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts3345[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21j = &(this->inputsTPParent->u21j_darts3345[this->getLocalID()]);
(void)u21j/*OMP_SHARED_PRIVATE*/;
double** u21jm1 = &(this->inputsTPParent->u21jm1_darts3345[this->getLocalID()]);
(void)u21jm1/*OMP_SHARED_PRIVATE*/;
double** u31j = &(this->inputsTPParent->u31j_darts3345[this->getLocalID()]);
(void)u31j/*OMP_SHARED_PRIVATE*/;
double** u31jm1 = &(this->inputsTPParent->u31jm1_darts3345[this->getLocalID()]);
(void)u31jm1/*OMP_SHARED_PRIVATE*/;
double** u41j = &(this->inputsTPParent->u41j_darts3345[this->getLocalID()]);
(void)u41j/*OMP_SHARED_PRIVATE*/;
double** u41jm1 = &(this->inputsTPParent->u41jm1_darts3345[this->getLocalID()]);
(void)u41jm1/*OMP_SHARED_PRIVATE*/;
double** u51j = &(this->inputsTPParent->u51j_darts3345[this->getLocalID()]);
(void)u51j/*OMP_SHARED_PRIVATE*/;
double** u51jm1 = &(this->inputsTPParent->u51jm1_darts3345[this->getLocalID()]);
(void)u51jm1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3345((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3345[0].decDep();
return;
}
for (int i_darts_counter_temp3345 = (*i);i_darts_counter_temp3345<=endRange && i_darts_counter_temp3345<=this->inputsTPParent->lastIteration3345;i_darts_counter_temp3345++)
{
{
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp3345 = (*k);
for(;k_darts_counter_temp3345 <= nz - 2;k_darts_counter_temp3345++){
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp3345 = (*j);
for(;j_darts_counter_temp3345 <= jend;j_darts_counter_temp3345++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3345 = (*m);
for(;m_darts_counter_temp3345 < 5;m_darts_counter_temp3345++){
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][m_darts_counter_temp3345] - ty2 * (flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][m_darts_counter_temp3345] - flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][m_darts_counter_temp3345]);
}
(*m) = m_darts_counter_temp3345;
}
}
(*j) = j_darts_counter_temp3345;
}
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp3345 = (*j);
for(;j_darts_counter_temp3345 <= (*(*L2));j_darts_counter_temp3345++){
(*(*tmp)) = 1. / rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][0];
(*(*u21j)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1];
(*(*u31j)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2];
(*(*u41j)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3];
(*(*u51j)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4];
(*(*tmp)) = 1. / rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][0];
(*(*u21jm1)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][1];
(*(*u31jm1)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][2];
(*(*u41jm1)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][3];
(*(*u51jm1)) = (*(*tmp)) * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][4];
flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1] = ty3 * ((*(*u21j)) - (*(*u21jm1)));
flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2] = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3] = ty3 * ((*(*u41j)) - (*(*u41jm1)));
flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * ty3 * (((*(*u21j)) * (*(*u21j)) + (*(*u31j)) * (*(*u31j)) + (*(*u41j)) * (*(*u41j))) - ((*(*u21jm1)) * (*(*u21jm1)) + (*(*u31jm1)) * (*(*u31jm1)) + (*(*u41jm1)) * (*(*u41jm1)))) + (1. / 6.) * ty3 * ((*(*u31j)) * (*(*u31j)) - (*(*u31jm1)) * (*(*u31jm1))) + 1.3999999999999999 * 1.3999999999999999 * ty3 * ((*(*u51j)) - (*(*u51jm1)));
}
(*j) = j_darts_counter_temp3345;
}
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp3345 = (*j);
for(;j_darts_counter_temp3345 <= jend;j_darts_counter_temp3345++){
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][0] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][0] + dy1 * ty1 * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][0] - 2. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][0] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][0]);
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][1] - flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1]) + dy2 * ty1 * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][1] - 2. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][1] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][1]);
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][2] - flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2]) + dy3 * ty1 * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][2] - 2. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][2] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][2]);
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][3] - flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3]) + dy4 * ty1 * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][3] - 2. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][3] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][3]);
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][4] - flux[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4]) + dy5 * ty1 * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][4] - 2. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][4] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][4]);
}
(*j) = j_darts_counter_temp3345;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3345 = (*m);
for(;m_darts_counter_temp3345 < 5;m_darts_counter_temp3345++){
frct[(i_darts_counter_temp3345)][1][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][1][k_darts_counter_temp3345][m_darts_counter_temp3345] - (*(*dsspm)) * (+5. * rsd[(i_darts_counter_temp3345)][1][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][2][k_darts_counter_temp3345][m_darts_counter_temp3345] + rsd[(i_darts_counter_temp3345)][3][k_darts_counter_temp3345][m_darts_counter_temp3345]);
frct[(i_darts_counter_temp3345)][2][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][2][k_darts_counter_temp3345][m_darts_counter_temp3345] - (*(*dsspm)) * (-4. * rsd[(i_darts_counter_temp3345)][1][k_darts_counter_temp3345][m_darts_counter_temp3345] + 6. * rsd[(i_darts_counter_temp3345)][2][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][3][k_darts_counter_temp3345][m_darts_counter_temp3345] + rsd[(i_darts_counter_temp3345)][4][k_darts_counter_temp3345][m_darts_counter_temp3345]);
}
(*m) = m_darts_counter_temp3345;
}
(*(*jst1)) = 3;
(*(*jend1)) = ny - 4;
{
/*Loop's init*/
(*j) = (*(*jst1));
int j_darts_counter_temp3345 = (*j);
for(;j_darts_counter_temp3345 <= (*(*jend1));j_darts_counter_temp3345++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3345 = (*m);
for(;m_darts_counter_temp3345 < 5;m_darts_counter_temp3345++){
frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][m_darts_counter_temp3345] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 2][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 - 1][k_darts_counter_temp3345][m_darts_counter_temp3345] + 6. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 1][k_darts_counter_temp3345][m_darts_counter_temp3345] + rsd[(i_darts_counter_temp3345)][j_darts_counter_temp3345 + 2][k_darts_counter_temp3345][m_darts_counter_temp3345]);
}
(*m) = m_darts_counter_temp3345;
}
}
(*j) = j_darts_counter_temp3345;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3345 = (*m);
for(;m_darts_counter_temp3345 < 5;m_darts_counter_temp3345++){
frct[(i_darts_counter_temp3345)][ny - 3][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][ny - 3][k_darts_counter_temp3345][m_darts_counter_temp3345] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3345)][ny - 5][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][ny - 4][k_darts_counter_temp3345][m_darts_counter_temp3345] + 6. * rsd[(i_darts_counter_temp3345)][ny - 3][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][ny - 2][k_darts_counter_temp3345][m_darts_counter_temp3345]);
frct[(i_darts_counter_temp3345)][ny - 2][k_darts_counter_temp3345][m_darts_counter_temp3345] = frct[(i_darts_counter_temp3345)][ny - 2][k_darts_counter_temp3345][m_darts_counter_temp3345] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3345)][ny - 4][k_darts_counter_temp3345][m_darts_counter_temp3345] - 4. * rsd[(i_darts_counter_temp3345)][ny - 3][k_darts_counter_temp3345][m_darts_counter_temp3345] + 5. * rsd[(i_darts_counter_temp3345)][ny - 2][k_darts_counter_temp3345][m_darts_counter_temp3345]);
}
(*m) = m_darts_counter_temp3345;
}
}
(*k) = k_darts_counter_temp3345;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3345[0].decDep();
}
TP3345::TP3345(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3345** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L2_darts3345(new int*[this->numThreads]),dsspm_darts3345(new double*[this->numThreads]),i_darts3345(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts3345(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jend1_darts3345(new int*[this->numThreads]),jst1_darts3345(new int*[this->numThreads]),k_darts3345(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts3345(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts3345(new double*[this->numThreads]),u21j_darts3345(new double*[this->numThreads]),u21jm1_darts3345(new double*[this->numThreads]),u31j_darts3345(new double*[this->numThreads]),u31jm1_darts3345(new double*[this->numThreads]),u41j_darts3345(new double*[this->numThreads]),u41jm1_darts3345(new double*[this->numThreads]),u51j_darts3345(new double*[this->numThreads]),u51jm1_darts3345(new double*[this->numThreads]), initIteration3345(in_initIteration), lastIteration3345(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets3345(new _barrierCodelets3345[1]) ,checkInCodelets3346(new _checkInCodelets3346[this->numThreads]){
/*Initialize the loop parameters*/
range3345 = abs (lastIteration3345 - initIteration3345) / 1;
rangePerCodelet3345 = range3345 / numThreads;
minIteration3345 = min<int>(lastIteration3345, initIteration3345);
remainderRange3345 = range3345 % numThreads;
/*Initialize inputs and vars.*/
this->L2_darts3345 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->dsspm_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jend1_darts3345 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jst1_darts3345 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21j_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21jm1_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31j_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31jm1_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41j_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41jm1_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51j_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51jm1_darts3345 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets3345[0] = _barrierCodelets3345(this->numThreads,this->numThreads,this, 0);
_checkInCodelets3346 * checkInCodelets3346Ptr = (this->checkInCodelets3346);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets3346);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets3346Ptr) = _checkInCodelets3346(2,1,this,codeletCounter);
#else
(*checkInCodelets3346Ptr) = _checkInCodelets3346(1,1,this,codeletCounter);
#endif
(*checkInCodelets3346Ptr).decDep();
checkInCodelets3346Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP3345::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets3346[localID].setID (codeletID);
this->checkInCodelets3346[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets3346[localID + this->baseNumThreads * i] = _checkInCodelets3346(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets3346[localID + this->baseNumThreads * i] = _checkInCodelets3346(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets3346[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets3346[localID + this->baseNumThreads * i].decDep();
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
TP3345::~TP3345(){
delete [] L2_darts3345;
delete [] dsspm_darts3345;
delete [] jend1_darts3345;
delete [] jst1_darts3345;
delete [] tmp_darts3345;
delete [] u21j_darts3345;
delete [] u21jm1_darts3345;
delete [] u31j_darts3345;
delete [] u31jm1_darts3345;
delete [] u41j_darts3345;
delete [] u41jm1_darts3345;
delete [] u51j_darts3345;
delete [] u51jm1_darts3345;
delete [] barrierCodelets3345;
delete [] checkInCodelets3346;
}
/*TP3980: OMPForDirective*/
void TP3980::_barrierCodelets3980::fire(void)
{
TP3980* myTP = static_cast<TP3980*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets3980[0].decDep ();
}
bool TP3980::requestNewRangeIterations3980(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet3980 * codeletID;
int tempEndRange   = rangePerCodelet3980 * (codeletID + 1);
if (remainderRange3980 != 0)
{
if (codeletID < (uint32_t)remainderRange3980)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange3980;
tempEndRange += remainderRange3980;
}
}
tempStartRange = tempStartRange*1 + minIteration3980;
tempEndRange = tempEndRange*1 + minIteration3980;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration3980 < lastIteration3980)
{
(this->inputsTPParent->i_darts3980[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts3980[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration3980;
}
}
return isThereNewIteration;
}
void TP3980::_checkInCodelets3981::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->dsspm_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21k_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21km1_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31k_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31km1_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41k_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41km1_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51k_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts2204[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51km1_darts3980[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts2204[this->getID()]);

/*printing node 3981: ForStmt*/
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
double** dsspm = &(this->inputsTPParent->dsspm_darts3980[this->getLocalID()]);
(void)dsspm/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts3980[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts3980[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts3980[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts3980[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts3980[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts3980[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21k = &(this->inputsTPParent->u21k_darts3980[this->getLocalID()]);
(void)u21k/*OMP_SHARED_PRIVATE*/;
double** u21km1 = &(this->inputsTPParent->u21km1_darts3980[this->getLocalID()]);
(void)u21km1/*OMP_SHARED_PRIVATE*/;
double** u31k = &(this->inputsTPParent->u31k_darts3980[this->getLocalID()]);
(void)u31k/*OMP_SHARED_PRIVATE*/;
double** u31km1 = &(this->inputsTPParent->u31km1_darts3980[this->getLocalID()]);
(void)u31km1/*OMP_SHARED_PRIVATE*/;
double** u41 = &(this->inputsTPParent->u41_darts3980[this->getLocalID()]);
(void)u41/*OMP_SHARED_PRIVATE*/;
double** u41k = &(this->inputsTPParent->u41k_darts3980[this->getLocalID()]);
(void)u41k/*OMP_SHARED_PRIVATE*/;
double** u41km1 = &(this->inputsTPParent->u41km1_darts3980[this->getLocalID()]);
(void)u41km1/*OMP_SHARED_PRIVATE*/;
double** u51k = &(this->inputsTPParent->u51k_darts3980[this->getLocalID()]);
(void)u51k/*OMP_SHARED_PRIVATE*/;
double** u51km1 = &(this->inputsTPParent->u51km1_darts3980[this->getLocalID()]);
(void)u51km1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3980((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3980[0].decDep();
return;
}
for (int i_darts_counter_temp3980 = (*i);i_darts_counter_temp3980<=endRange && i_darts_counter_temp3980<=this->inputsTPParent->lastIteration3980;i_darts_counter_temp3980++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp3980 = (*j);
for(;j_darts_counter_temp3980 <= jend;j_darts_counter_temp3980++){
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp3980 = (*k);
for(;k_darts_counter_temp3980 <= nz - 1;k_darts_counter_temp3980++){
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0] = rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3];
(*(*u41)) = rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] / rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0];
(*(*q)) = 0.5 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] + rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] + rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3]) / rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0];
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] = rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] * (*(*u41));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] = rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] * (*(*u41));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] = rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] * (*(*u41)) + 0.40000000000000002 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] - (*(*q)));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] = (1.3999999999999999 * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] - 0.40000000000000002 * (*(*q))) * (*(*u41));
}
(*k) = k_darts_counter_temp3980;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp3980 = (*k);
for(;k_darts_counter_temp3980 <= nz - 2;k_darts_counter_temp3980++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3980 = (*m);
for(;m_darts_counter_temp3980 < 5;m_darts_counter_temp3980++){
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][m_darts_counter_temp3980] - tz2 * (flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][m_darts_counter_temp3980] - flux[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][m_darts_counter_temp3980]);
}
(*m) = m_darts_counter_temp3980;
}
}
(*k) = k_darts_counter_temp3980;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp3980 = (*k);
for(;k_darts_counter_temp3980 <= nz - 1;k_darts_counter_temp3980++){
(*(*tmp)) = 1. / rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0];
(*(*u21k)) = (*(*tmp)) * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1];
(*(*u31k)) = (*(*tmp)) * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2];
(*(*u41k)) = (*(*tmp)) * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3];
(*(*u51k)) = (*(*tmp)) * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4];
(*(*tmp)) = 1. / rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][0];
(*(*u21km1)) = (*(*tmp)) * rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][1];
(*(*u31km1)) = (*(*tmp)) * rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][2];
(*(*u41km1)) = (*(*tmp)) * rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][3];
(*(*u51km1)) = (*(*tmp)) * rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][4];
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] = tz3 * ((*(*u21k)) - (*(*u21km1)));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] = tz3 * ((*(*u31k)) - (*(*u31km1)));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tz3 * (((*(*u21k)) * (*(*u21k)) + (*(*u31k)) * (*(*u31k)) + (*(*u41k)) * (*(*u41k))) - ((*(*u21km1)) * (*(*u21km1)) + (*(*u31km1)) * (*(*u31km1)) + (*(*u41km1)) * (*(*u41km1)))) + (1. / 6.) * tz3 * ((*(*u41k)) * (*(*u41k)) - (*(*u41km1)) * (*(*u41km1))) + 1.3999999999999999 * 1.3999999999999999 * tz3 * ((*(*u51k)) - (*(*u51km1)));
}
(*k) = k_darts_counter_temp3980;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp3980 = (*k);
for(;k_darts_counter_temp3980 <= nz - 2;k_darts_counter_temp3980++){
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0] + dz1 * tz1 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][0] - 2. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][0] + rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][0]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][1] - flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1]) + dz2 * tz1 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][1] - 2. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][1] + rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][1]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][2] - flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2]) + dz3 * tz1 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][2] - 2. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][2] + rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][2]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][3] - flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3]) + dz4 * tz1 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][3] - 2. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][3] + rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][3]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][4] - flux[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4]) + dz5 * tz1 * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][4] - 2. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][4] + rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][4]);
}
(*k) = k_darts_counter_temp3980;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3980 = (*m);
for(;m_darts_counter_temp3980 < 5;m_darts_counter_temp3980++){
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][1][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][1][m_darts_counter_temp3980] - (*(*dsspm)) * (+5. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][1][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][2][m_darts_counter_temp3980] + rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][3][m_darts_counter_temp3980]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][2][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][2][m_darts_counter_temp3980] - (*(*dsspm)) * (-4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][1][m_darts_counter_temp3980] + 6. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][2][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][3][m_darts_counter_temp3980] + rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][4][m_darts_counter_temp3980]);
}
(*m) = m_darts_counter_temp3980;
}
{
/*Loop's init*/
(*k) = 3;
int k_darts_counter_temp3980 = (*k);
for(;k_darts_counter_temp3980 <= nz - 4;k_darts_counter_temp3980++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3980 = (*m);
for(;m_darts_counter_temp3980 < 5;m_darts_counter_temp3980++){
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][m_darts_counter_temp3980] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 - 2][m_darts_counter_temp3980] - 4. * rsd[k_darts_counter_temp3980 - 1][(i_darts_counter_temp3980)][j_darts_counter_temp3980][m_darts_counter_temp3980] + 6. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 1][m_darts_counter_temp3980] + rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][k_darts_counter_temp3980 + 2][m_darts_counter_temp3980]);
}
(*m) = m_darts_counter_temp3980;
}
}
(*k) = k_darts_counter_temp3980;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp3980 = (*m);
for(;m_darts_counter_temp3980 < 5;m_darts_counter_temp3980++){
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 3][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 3][m_darts_counter_temp3980] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 5][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 4][m_darts_counter_temp3980] + 6. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 3][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 2][m_darts_counter_temp3980]);
frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 2][m_darts_counter_temp3980] = frct[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 2][m_darts_counter_temp3980] - (*(*dsspm)) * (rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 4][m_darts_counter_temp3980] - 4. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 3][m_darts_counter_temp3980] + 5. * rsd[(i_darts_counter_temp3980)][j_darts_counter_temp3980][nz - 2][m_darts_counter_temp3980]);
}
(*m) = m_darts_counter_temp3980;
}
}
(*j) = j_darts_counter_temp3980;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets3980[0].decDep();
}
TP3980::TP3980(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3980** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),dsspm_darts3980(new double*[this->numThreads]),i_darts3980(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts3980(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts3980(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts3980(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts3980(new double*[this->numThreads]),tmp_darts3980(new double*[this->numThreads]),u21k_darts3980(new double*[this->numThreads]),u21km1_darts3980(new double*[this->numThreads]),u31k_darts3980(new double*[this->numThreads]),u31km1_darts3980(new double*[this->numThreads]),u41_darts3980(new double*[this->numThreads]),u41k_darts3980(new double*[this->numThreads]),u41km1_darts3980(new double*[this->numThreads]),u51k_darts3980(new double*[this->numThreads]),u51km1_darts3980(new double*[this->numThreads]), initIteration3980(in_initIteration), lastIteration3980(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets3980(new _barrierCodelets3980[1]) ,checkInCodelets3981(new _checkInCodelets3981[this->numThreads]){
/*Initialize the loop parameters*/
range3980 = abs (lastIteration3980 - initIteration3980) / 1;
rangePerCodelet3980 = range3980 / numThreads;
minIteration3980 = min<int>(lastIteration3980, initIteration3980);
remainderRange3980 = range3980 % numThreads;
/*Initialize inputs and vars.*/
this->dsspm_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->q_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21k_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21km1_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31k_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31km1_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41k_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41km1_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51k_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51km1_darts3980 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets3980[0] = _barrierCodelets3980(this->numThreads,this->numThreads,this, 0);
_checkInCodelets3981 * checkInCodelets3981Ptr = (this->checkInCodelets3981);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets3981);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets3981Ptr) = _checkInCodelets3981(2,1,this,codeletCounter);
#else
(*checkInCodelets3981Ptr) = _checkInCodelets3981(1,1,this,codeletCounter);
#endif
(*checkInCodelets3981Ptr).decDep();
checkInCodelets3981Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP3980::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets3981[localID].setID (codeletID);
this->checkInCodelets3981[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets3981[localID + this->baseNumThreads * i] = _checkInCodelets3981(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets3981[localID + this->baseNumThreads * i] = _checkInCodelets3981(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets3981[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets3981[localID + this->baseNumThreads * i].decDep();
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
TP3980::~TP3980(){
delete [] dsspm_darts3980;
delete [] q_darts3980;
delete [] tmp_darts3980;
delete [] u21k_darts3980;
delete [] u21km1_darts3980;
delete [] u31k_darts3980;
delete [] u31km1_darts3980;
delete [] u41_darts3980;
delete [] u41k_darts3980;
delete [] u41km1_darts3980;
delete [] u51k_darts3980;
delete [] u51km1_darts3980;
delete [] barrierCodelets3980;
delete [] checkInCodelets3981;
}
/*TP7: TP_jacld*/
void TP7::_checkInCodelets4885::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 4885: DeclStmt*/

/*printing node 4886: DeclStmt*/

/*printing node 4887: DeclStmt*/

/*printing node 4888: DeclStmt*/

/*printing node 4889: DeclStmt*/

/*printing node 4890: BinaryOperator*/
(this->inputsTPParent->r43_darts7[this->getID()]) = (4. / 3.);

/*printing node 4894: BinaryOperator*/
(this->inputsTPParent->c1345_darts7[this->getID()]) = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

/*printing node 4902: BinaryOperator*/
(this->inputsTPParent->c34_darts7[this->getID()]) = 0.10000000000000001 * 1.;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 4885 nextRegion: 4906 */
myTP->controlTPParent->checkInCodelets4906[this->getID()].decDep();
}
void TP7::_checkInCodelets4906::fire(void)
{
/*region 4906 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP4906;
if(idx < myTP->TPsToUse4906){
if (!__sync_val_compare_and_swap (&(myTP->TP4906_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse4906;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse4906;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse4906 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse4906 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP4906 > (myTP, myTP->codeletsPerTP4906 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP4906Ptr[idx]));
#else
place < TP4906 > (idx, myTP, myTP->codeletsPerTP4906 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP4906Ptr[idx]));
#endif
}else{
if (myTP->TP4906Ptr[idx] != nullptr){
myTP->TP4906Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP7::_barrierCodelets4906::fire(void)
{
TP7* myTP =  static_cast<TP7*>(myTP_);

for(size_t codeletsCounter=0; codeletsCounter < (size_t)myTP->numThreads; codeletsCounter++)
{
myTP->nextCodeletsjacld[codeletsCounter]->decDep();
}
}
TP7::TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k):ompTP(in_numThreads, in_mainCodeletID), ptrToThisFunctionTP(in_ptrToThisFunctionTP), inputsTPParent(this), controlTPParent(this), nextCodeletsjacld( new Codelet*[in_numThreads])
, nextSyncCodeletsjacld( new Codelet*[in_numThreads])
,k_darts7(new int[this->numThreads])
,c1345_darts7(new double[this->numThreads])
,c34_darts7(new double[this->numThreads])
,i_darts7(new int[this->numThreads])
,j_darts7(new int[this->numThreads])
,r43_darts7(new double[this->numThreads])
,tmp1_darts7(new double[this->numThreads])
,tmp2_darts7(new double[this->numThreads])
,tmp3_darts7(new double[this->numThreads])
, TP4906Ptr(new TP4906 *[NUMTPS4906]), TP4906_alreadyLaunched(new size_t [NUMTPS4906]), numTPsSet4906(0), numTPsReady4906(0), TPsToUse4906(NUMTPS4906), codeletsPerTP4906(this->numThreads/NUMTPS4906), totalCodelets4906(this->TPsToUse4906*this->codeletsPerTP4906) ,checkInCodelets4885(new _checkInCodelets4885[this->numThreads]) ,checkInCodelets4906(new _checkInCodelets4906[this->numThreads]) ,barrierCodelets4906(new _barrierCodelets4906[1]){
barrierCodelets4906[0] = _barrierCodelets4906(NUMTPS4906,NUMTPS4906,this, 0);
_checkInCodelets4906 * checkInCodelets4906Ptr = (this->checkInCodelets4906);
for(int i=0; i<NUMTPS4906; i++)
{
TP4906Ptr[i] = nullptr;
TP4906_alreadyLaunched[i] = 0;
}
_checkInCodelets4885 * checkInCodelets4885Ptr = (this->checkInCodelets4885);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets4885);
#endif
for(size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++)
{
(*checkInCodelets4906Ptr) = _checkInCodelets4906(1,1,this,codeletCounter);
checkInCodelets4906Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets4885Ptr) = _checkInCodelets4885(2,1,this,codeletCounter);
#else
(*checkInCodelets4885Ptr) = _checkInCodelets4885(1,1,this,codeletCounter);
#endif
(*checkInCodelets4885Ptr).decDep();
checkInCodelets4885Ptr++;
}
if(this->numThreads == 1){
this->nextCodeletsjacld[0] = in_mainNextCodelet;
this->nextSyncCodeletsjacld[0] = in_mainSyncCodelet;
this->k_darts7[0]= in_k;
this->availableCodelets[0] = 1;
}
else
{
this->k_darts7[this->mainCodeletID]= in_k;
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
TP7::~TP7(){
delete [] barrierCodelets4906;
delete [] checkInCodelets4906;
delete [] checkInCodelets4885;
delete [] nextCodeletsjacld;
delete [] nextSyncCodeletsjacld;
delete [] k_darts7;
delete [] c1345_darts7;
delete [] c34_darts7;
delete [] i_darts7;
delete [] j_darts7;
delete [] r43_darts7;
delete [] tmp1_darts7;
delete [] tmp2_darts7;
delete [] tmp3_darts7;
}
void TP7::setNewInputs(int in_k, size_t codeletID){
this->k_darts7[codeletID]= in_k;
}
/*TP4906: OMPForDirective*/
void TP4906::_barrierCodelets4906::fire(void)
{
TP4906* myTP = static_cast<TP4906*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets4906[0].decDep ();
}
bool TP4906::requestNewRangeIterations4906(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 1*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet4906 * codeletID;
int tempEndRange   = rangePerCodelet4906 * (codeletID + 1);
if (remainderRange4906 != 0)
{
if (codeletID < (uint32_t)remainderRange4906)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange4906;
tempEndRange += remainderRange4906;
}
}
tempStartRange = tempStartRange*1 + minIteration4906;
tempEndRange = tempEndRange*1 + minIteration4906;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration4906 < lastIteration4906)
{
(this->inputsTPParent->i_darts4906[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts4906[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration4906;
}
}
return isThereNewIteration;
}
void TP4906::_checkInCodelets4907::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->c1345_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->c34_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->c34_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->k_darts4906[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->k_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->r43_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->r43_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp1_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp2_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts7[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp3_darts4906[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts7[this->getID()]);

/*printing node 4907: ForStmt*/
/*var: c1345*/
/*var: c34*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: r43*/
/*var: tmp1*/
/*var: tmp2*/
/*var: tmp3*/
double** c1345 = &(this->inputsTPParent->c1345_darts4906[this->getLocalID()]);
(void)c1345/*OMP_SHARED_PRIVATE*/;
double** c34 = &(this->inputsTPParent->c34_darts4906[this->getLocalID()]);
(void)c34/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts4906[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts4906[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** k = &(this->inputsTPParent->k_darts4906[this->getLocalID()]);
(void)k/*OMP_SHARED_PRIVATE*/;
double** r43 = &(this->inputsTPParent->r43_darts4906[this->getLocalID()]);
(void)r43/*OMP_SHARED_PRIVATE*/;
double** tmp1 = &(this->inputsTPParent->tmp1_darts4906[this->getLocalID()]);
(void)tmp1/*OMP_SHARED_PRIVATE*/;
double** tmp2 = &(this->inputsTPParent->tmp2_darts4906[this->getLocalID()]);
(void)tmp2/*OMP_SHARED_PRIVATE*/;
double** tmp3 = &(this->inputsTPParent->tmp3_darts4906[this->getLocalID()]);
(void)tmp3/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations4906((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets4906[0].decDep();
return;
}
for (int i_darts_counter_temp4906 = (*i);i_darts_counter_temp4906<=endRange && i_darts_counter_temp4906<=this->inputsTPParent->lastIteration4906;i_darts_counter_temp4906++)
{
{
int id  = omp_get_thread_num();
printf("Thread #%d: jacld entering loopi %d\n", id, (i_darts_counter_temp4906));
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp4906 = (*j);
for(;j_darts_counter_temp4906 <= jend;j_darts_counter_temp4906++){
(*(*tmp1)) = 1. / u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][0] = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][1] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][2] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][3] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][4] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][0] = dt * 2. * (tx1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]) + ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]) + tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]));
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][1] = 1. + dt * 2. * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][2] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][3] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][4] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][0] = dt * 2. * (tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]) + ty1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]) + tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]));
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][1] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][2] = 1. + dt * 2. * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][3] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][4] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][0] = dt * 2. * (tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]) + ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]) + tz1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]));
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][1] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][2] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][3] = 1. + dt * 2. * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][4] = 0.;
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][0] = dt * 2. * (tx1 * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][4]) + ty1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][4]) + tz1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][4]));
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][1] = dt * 2. * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1] + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1] + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][1]);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][2] = dt * 2. * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2] + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2] + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][2]);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][3] = dt * 2. * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3] + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3] + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906][(*(*k))][3]);
d[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][4] = 1. + dt * 2. * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1)) + tz1 * (*(*c1345)) * (*(*tmp1))) + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
(*(*tmp1)) = 1. / u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][0] = -dt * tz1 * dz1;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][1] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][2] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][3] = -dt * tz2;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][4] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][0] = -dt * tz2 * (-(u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2))) - dt * tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1]);
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][1] = -dt * tz2 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1))) - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][2] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][3] = -dt * tz2 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * (*(*tmp1)));
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][4] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][0] = -dt * tz2 * (-(u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2))) - dt * tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2]);
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][1] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][2] = -dt * tz2 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1))) - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][3] = -dt * tz2 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * (*(*tmp1)));
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][4] = 0.;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][0] = -dt * tz2 * (-(u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1))) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1))) + 0.5 * 0.40000000000000002 * ((u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] + u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] + u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2)))) - dt * tz1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]);
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][1] = -dt * tz2 * (-0.40000000000000002 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * (*(*tmp1))));
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][2] = -dt * tz2 * (-0.40000000000000002 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * (*(*tmp1))));
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][3] = -dt * tz2 * (2. - 0.40000000000000002) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1))) - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][4] = -dt * tz2 * 0.40000000000000002;
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][0] = -dt * tz2 * ((0.40000000000000002 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] + u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] + u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][4] * (*(*tmp1)))) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1)))) - dt * tz1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1]) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2]) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) - (*(*c1345)) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][4]);
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][1] = -dt * tz2 * (-0.40000000000000002 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2))) - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1];
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][2] = -dt * tz2 * (-0.40000000000000002 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2))) - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2];
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][3] = -dt * tz2 * (1.3999999999999999 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][1] + u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][2] + 3. * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3]) * (*(*tmp2)))) - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3];
a[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][4] = -dt * tz2 * (1.3999999999999999 * (u[(*(*k)) - 1][(i_darts_counter_temp4906)][j_darts_counter_temp4906][3] * (*(*tmp1)))) - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
(*(*tmp1)) = 1. / u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][0] = -dt * ty1 * dy1;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][1] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][2] = -dt * ty2;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][3] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][4] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][0] = -dt * ty2 * (-(u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2]) * (*(*tmp2))) - dt * ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1]);
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][1] = -dt * ty2 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1))) - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][2] = -dt * ty2 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * (*(*tmp1)));
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][3] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][4] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][0] = -dt * ty2 * (-(u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1))) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1))) + 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] + u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] + u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (*(*tmp2)))) - dt * ty1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2]);
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][1] = -dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * (*(*tmp1))));
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][2] = -dt * ty2 * ((2. - 0.40000000000000002) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][3] = -dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3] * (*(*tmp1))));
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][4] = -dt * ty2 * 0.40000000000000002;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][0] = -dt * ty2 * (-(u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (*(*tmp2))) - dt * ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]);
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][1] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][2] = -dt * ty2 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3] * (*(*tmp1)));
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][3] = -dt * ty2 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1))) - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][4] = 0.;
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][0] = -dt * ty2 * ((0.40000000000000002 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] + u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] + u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][4] * (*(*tmp1)))) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]))) - (*(*c1345)) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][4]);
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][1] = -dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2]) * (*(*tmp2))) - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1];
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][2] = -dt * ty2 * (1.3999999999999999 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][1] + 3. * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] + u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (*(*tmp2)))) - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2];
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][3] = -dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3]) * (*(*tmp2))) - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][3];
b[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][4] = -dt * ty2 * (1.3999999999999999 * (u[(i_darts_counter_temp4906)][j_darts_counter_temp4906 - 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
(*(*tmp1)) = 1. / u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][0] = -dt * tx1 * dx1;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][1] = -dt * tx2;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][2] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][3] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][0][4] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][0] = -dt * tx2 * (-(u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1))) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1))) + 0.40000000000000002 * 0.5 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]) * (*(*tmp2))) - dt * tx1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1]);
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][1] = -dt * tx2 * ((2. - 0.40000000000000002) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][2] = -dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * (*(*tmp1))));
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][3] = -dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * (*(*tmp1))));
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][1][4] = -dt * tx2 * 0.40000000000000002;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][0] = -dt * tx2 * (-(u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2]) * (*(*tmp2))) - dt * tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2]);
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][1] = -dt * tx2 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * (*(*tmp1)));
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][2] = -dt * tx2 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1))) - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][3] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][2][4] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][0] = -dt * tx2 * (-(u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]) * (*(*tmp2))) - dt * tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]);
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][1] = -dt * tx2 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * (*(*tmp1)));
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][2] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][3] = -dt * tx2 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1))) - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][3][4] = 0.;
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][0] = -dt * tx2 * ((0.40000000000000002 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][4] * (*(*tmp1)))) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1]) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2]) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]) * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]))) - (*(*c1345)) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][4]);
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][1] = -dt * tx2 * (1.3999999999999999 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((3. * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] + u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3]) * (*(*tmp2)))) - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1];
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][2] = -dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1]) * (*(*tmp2))) - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][2];
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][3] = -dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3] * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1]) * (*(*tmp2))) - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][3];
c[(i_darts_counter_temp4906)][j_darts_counter_temp4906][4][4] = -dt * tx2 * (1.3999999999999999 * (u[(i_darts_counter_temp4906) - 1][j_darts_counter_temp4906][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
}
(*j) = j_darts_counter_temp4906;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets4906[0].decDep();
}
TP4906::TP4906(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration, int in_lastIteration, TP4906** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),c1345_darts4906(new double*[this->numThreads]),c34_darts4906(new double*[this->numThreads]),i_darts4906(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts4906(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts4906(new int*[this->numThreads]),r43_darts4906(new double*[this->numThreads]),tmp1_darts4906(new double*[this->numThreads]),tmp2_darts4906(new double*[this->numThreads]),tmp3_darts4906(new double*[this->numThreads]), initIteration4906(in_initIteration), lastIteration4906(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets4906(new _barrierCodelets4906[1]) ,checkInCodelets4907(new _checkInCodelets4907[this->numThreads]){
/*Initialize the loop parameters*/
range4906 = abs (lastIteration4906 - initIteration4906) / 1;
rangePerCodelet4906 = range4906 / numThreads;
minIteration4906 = min<int>(lastIteration4906, initIteration4906);
remainderRange4906 = range4906 % numThreads;
/*Initialize inputs and vars.*/
this->c1345_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->c34_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->k_darts4906 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->r43_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp1_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp2_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp3_darts4906 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets4906[0] = _barrierCodelets4906(this->numThreads,this->numThreads,this, 0);
_checkInCodelets4907 * checkInCodelets4907Ptr = (this->checkInCodelets4907);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets4907);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets4907Ptr) = _checkInCodelets4907(2,1,this,codeletCounter);
#else
(*checkInCodelets4907Ptr) = _checkInCodelets4907(1,1,this,codeletCounter);
#endif
(*checkInCodelets4907Ptr).decDep();
checkInCodelets4907Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP4906::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets4907[localID].setID (codeletID);
this->checkInCodelets4907[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets4907[localID + this->baseNumThreads * i] = _checkInCodelets4907(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets4907[localID + this->baseNumThreads * i] = _checkInCodelets4907(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets4907[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets4907[localID + this->baseNumThreads * i].decDep();
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
TP4906::~TP4906(){
delete [] c1345_darts4906;
delete [] c34_darts4906;
delete [] k_darts4906;
delete [] r43_darts4906;
delete [] tmp1_darts4906;
delete [] tmp2_darts4906;
delete [] tmp3_darts4906;
delete [] barrierCodelets4906;
delete [] checkInCodelets4907;
}
/*TP8: TP_jacu*/
void TP8::_checkInCodelets7407::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 7407: DeclStmt*/

/*printing node 7408: DeclStmt*/

/*printing node 7409: DeclStmt*/

/*printing node 7410: DeclStmt*/

/*printing node 7411: DeclStmt*/

/*printing node 7412: BinaryOperator*/
(this->inputsTPParent->r43_darts8[this->getID()]) = (4. / 3.);

/*printing node 7416: BinaryOperator*/
(this->inputsTPParent->c1345_darts8[this->getID()]) = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

/*printing node 7424: BinaryOperator*/
(this->inputsTPParent->c34_darts8[this->getID()]) = 0.10000000000000001 * 1.;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 7407 nextRegion: 7428 */
myTP->controlTPParent->checkInCodelets7428[this->getID()].decDep();
}
void TP8::_checkInCodelets7428::fire(void)
{
/*region 7428 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP7428;
if(idx < myTP->TPsToUse7428){
if (!__sync_val_compare_and_swap (&(myTP->TP7428_alreadyLaunched[idx]), 0, 1)){
int range = abs (ist - iend) / 1;
int rangePerCodelet = range / myTP->TPsToUse7428;
int minIteration = min<int >(ist, iend);
int remainderRange = range % myTP->TPsToUse7428;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(iend < ist)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == 0)
{
lastIteration = lastIteration - 1;
}
if(idx == myTP->TPsToUse7428 - 1)
{
lastIteration = ist;
}
#if USEINVOKE == 1
invoke < TP7428 > (myTP, myTP->codeletsPerTP7428 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP7428Ptr[idx]));
#else
place < TP7428 > (idx, myTP, myTP->codeletsPerTP7428 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP7428Ptr[idx]));
#endif
}else{
if (myTP->TP7428Ptr[idx] != nullptr){
myTP->TP7428Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP8::_barrierCodelets7428::fire(void)
{
TP8* myTP =  static_cast<TP8*>(myTP_);

for(size_t codeletsCounter=0; codeletsCounter < (size_t)myTP->numThreads; codeletsCounter++)
{
myTP->nextCodeletsjacu[codeletsCounter]->decDep();
}
}
TP8::TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k):ompTP(in_numThreads, in_mainCodeletID), ptrToThisFunctionTP(in_ptrToThisFunctionTP), inputsTPParent(this), controlTPParent(this), nextCodeletsjacu( new Codelet*[in_numThreads])
, nextSyncCodeletsjacu( new Codelet*[in_numThreads])
,k_darts8(new int[this->numThreads])
,c1345_darts8(new double[this->numThreads])
,c34_darts8(new double[this->numThreads])
,i_darts8(new int[this->numThreads])
,j_darts8(new int[this->numThreads])
,r43_darts8(new double[this->numThreads])
,tmp1_darts8(new double[this->numThreads])
,tmp2_darts8(new double[this->numThreads])
,tmp3_darts8(new double[this->numThreads])
, TP7428Ptr(new TP7428 *[NUMTPS7428]), TP7428_alreadyLaunched(new size_t [NUMTPS7428]), numTPsSet7428(0), numTPsReady7428(0), TPsToUse7428(NUMTPS7428), codeletsPerTP7428(this->numThreads/NUMTPS7428), totalCodelets7428(this->TPsToUse7428*this->codeletsPerTP7428) ,checkInCodelets7407(new _checkInCodelets7407[this->numThreads]) ,checkInCodelets7428(new _checkInCodelets7428[this->numThreads]) ,barrierCodelets7428(new _barrierCodelets7428[1]){
barrierCodelets7428[0] = _barrierCodelets7428(NUMTPS7428,NUMTPS7428,this, 0);
_checkInCodelets7428 * checkInCodelets7428Ptr = (this->checkInCodelets7428);
for(int i=0; i<NUMTPS7428; i++)
{
TP7428Ptr[i] = nullptr;
TP7428_alreadyLaunched[i] = 0;
}
_checkInCodelets7407 * checkInCodelets7407Ptr = (this->checkInCodelets7407);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets7407);
#endif
for(size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++)
{
(*checkInCodelets7428Ptr) = _checkInCodelets7428(1,1,this,codeletCounter);
checkInCodelets7428Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets7407Ptr) = _checkInCodelets7407(2,1,this,codeletCounter);
#else
(*checkInCodelets7407Ptr) = _checkInCodelets7407(1,1,this,codeletCounter);
#endif
(*checkInCodelets7407Ptr).decDep();
checkInCodelets7407Ptr++;
}
if(this->numThreads == 1){
this->nextCodeletsjacu[0] = in_mainNextCodelet;
this->nextSyncCodeletsjacu[0] = in_mainSyncCodelet;
this->k_darts8[0]= in_k;
this->availableCodelets[0] = 1;
}
else
{
this->k_darts8[this->mainCodeletID]= in_k;
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
TP8::~TP8(){
delete [] barrierCodelets7428;
delete [] checkInCodelets7428;
delete [] checkInCodelets7407;
delete [] nextCodeletsjacu;
delete [] nextSyncCodeletsjacu;
delete [] k_darts8;
delete [] c1345_darts8;
delete [] c34_darts8;
delete [] i_darts8;
delete [] j_darts8;
delete [] r43_darts8;
delete [] tmp1_darts8;
delete [] tmp2_darts8;
delete [] tmp3_darts8;
}
void TP8::setNewInputs(int in_k, size_t codeletID){
this->k_darts8[codeletID]= in_k;
}
/*TP7428: OMPForDirective*/
void TP7428::_barrierCodelets7428::fire(void)
{
TP7428* myTP = static_cast<TP7428*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets7428[0].decDep ();
}
bool TP7428::requestNewRangeIterations7428(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 1*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet7428 * codeletID;
int tempEndRange   = rangePerCodelet7428 * (codeletID + 1);
if (remainderRange7428 != 0)
{
if (codeletID < (uint32_t)remainderRange7428)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange7428;
tempEndRange += remainderRange7428;
}
}
tempStartRange = tempStartRange*1 + minIteration7428;
tempEndRange = tempEndRange*1 + minIteration7428;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration7428 < lastIteration7428)
{
(this->inputsTPParent->i_darts7428[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts7428[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == 0)
{
*endRange = *endRange - 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration7428;
}
}
return isThereNewIteration;
}
void TP7428::_checkInCodelets7429::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->c1345_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->c34_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->c34_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->k_darts7428[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->k_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->r43_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->r43_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp1_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp2_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts8[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp3_darts7428[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts8[this->getID()]);

/*printing node 7429: ForStmt*/
/*var: c1345*/
/*var: c34*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: r43*/
/*var: tmp1*/
/*var: tmp2*/
/*var: tmp3*/
double** c1345 = &(this->inputsTPParent->c1345_darts7428[this->getLocalID()]);
(void)c1345/*OMP_SHARED_PRIVATE*/;
double** c34 = &(this->inputsTPParent->c34_darts7428[this->getLocalID()]);
(void)c34/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts7428[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts7428[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** k = &(this->inputsTPParent->k_darts7428[this->getLocalID()]);
(void)k/*OMP_SHARED_PRIVATE*/;
double** r43 = &(this->inputsTPParent->r43_darts7428[this->getLocalID()]);
(void)r43/*OMP_SHARED_PRIVATE*/;
double** tmp1 = &(this->inputsTPParent->tmp1_darts7428[this->getLocalID()]);
(void)tmp1/*OMP_SHARED_PRIVATE*/;
double** tmp2 = &(this->inputsTPParent->tmp2_darts7428[this->getLocalID()]);
(void)tmp2/*OMP_SHARED_PRIVATE*/;
double** tmp3 = &(this->inputsTPParent->tmp3_darts7428[this->getLocalID()]);
(void)tmp3/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations7428((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets7428[0].decDep();
return;
}
for (int i_darts_counter_temp7428 = (*i);i_darts_counter_temp7428>=endRange && i_darts_counter_temp7428>=this->inputsTPParent->lastIteration7428;i_darts_counter_temp7428--)
{
{
{
/*Loop's init*/
(*j) = jend;
int j_darts_counter_temp7428 = (*j);
for(;j_darts_counter_temp7428 >= jst;j_darts_counter_temp7428--){
(*(*tmp1)) = 1. / u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][0] = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][1] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][2] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][3] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][4] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][0] = dt * 2. * (tx1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]) + ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]) + tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]));
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][1] = 1. + dt * 2. * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][2] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][3] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][4] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][0] = dt * 2. * (tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]) + ty1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]) + tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]));
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][1] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][2] = 1. + dt * 2. * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][3] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][4] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][0] = dt * 2. * (tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]) + ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]) + tz1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]));
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][1] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][2] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][3] = 1. + dt * 2. * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))) + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][4] = 0.;
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][0] = dt * 2. * (tx1 * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][4]) + ty1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][4]) + tz1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]))) - ((*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][4]));
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][1] = dt * 2. * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1] + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1] + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][1]);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][2] = dt * 2. * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2] + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2] + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][2]);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][3] = dt * 2. * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3] + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3] + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k))][3]);
d[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][4] = 1. + dt * 2. * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1)) + tz1 * (*(*c1345)) * (*(*tmp1))) + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
(*(*tmp1)) = 1. / u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][0] = -dt * tx1 * dx1;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][1] = dt * tx2;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][2] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][3] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][4] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][0] = dt * tx2 * (-(u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1))) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1))) + 0.40000000000000002 * 0.5 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]) * (*(*tmp2))) - dt * tx1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1]);
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][1] = dt * tx2 * ((2. - 0.40000000000000002) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][2] = dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * (*(*tmp1))));
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][3] = dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * (*(*tmp1))));
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][4] = dt * tx2 * 0.40000000000000002;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][0] = dt * tx2 * (-(u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2]) * (*(*tmp2))) - dt * tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2]);
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][1] = dt * tx2 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * (*(*tmp1)));
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][2] = dt * tx2 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1))) - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][3] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][4] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][0] = dt * tx2 * (-(u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]) * (*(*tmp2))) - dt * tx1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]);
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][1] = dt * tx2 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * (*(*tmp1)));
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][2] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][3] = dt * tx2 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1))) - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][4] = 0.;
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][0] = dt * tx2 * ((0.40000000000000002 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][4] * (*(*tmp1)))) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1]) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2]) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]) * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]))) - (*(*c1345)) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][4]);
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][1] = dt * tx2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((3. * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] + u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3]) * (*(*tmp2)))) - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1];
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][2] = dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1]) * (*(*tmp2))) - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][2];
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][3] = dt * tx2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3] * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1]) * (*(*tmp2))) - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][3];
a[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][4] = dt * tx2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428) + 1][j_darts_counter_temp7428][(*(*k))][1] * (*(*tmp1)))) - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
(*(*tmp1)) = 1. / u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][0] = -dt * ty1 * dy1;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][1] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][2] = dt * ty2;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][3] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][4] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][0] = dt * ty2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2]) * (*(*tmp2))) - dt * ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1]);
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][1] = dt * ty2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1))) - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][2] = dt * ty2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * (*(*tmp1)));
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][3] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][4] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][0] = dt * ty2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1))) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1))) + 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (*(*tmp2)))) - dt * ty1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2]);
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][1] = dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * (*(*tmp1))));
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][2] = dt * ty2 * ((2. - 0.40000000000000002) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][3] = dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3] * (*(*tmp1))));
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][4] = dt * ty2 * 0.40000000000000002;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][0] = dt * ty2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (*(*tmp2))) - dt * ty1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]);
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][1] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][2] = dt * ty2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3] * (*(*tmp1)));
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][3] = dt * ty2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1))) - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][4] = 0.;
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][0] = dt * ty2 * ((0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][4] * (*(*tmp1)))) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]))) - (*(*c1345)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][4]);
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][1] = dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2]) * (*(*tmp2))) - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1];
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][2] = dt * ty2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][1] + 3. * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (*(*tmp2)))) - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2];
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][3] = dt * ty2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3]) * (*(*tmp2))) - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][3];
b[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][4] = dt * ty2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428 + 1][(*(*k))][2] * (*(*tmp1)))) - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
(*(*tmp1)) = 1. / u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][0];
(*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
(*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][0] = -dt * tz1 * dz1;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][1] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][2] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][3] = dt * tz2;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][0][4] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][0] = dt * tz2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2))) - dt * tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1]);
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][1] = dt * tz2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1))) - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][2] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][3] = dt * tz2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * (*(*tmp1)));
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][1][4] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][0] = dt * tz2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2))) - dt * tz1 * (-(*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2]);
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][1] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][2] = dt * tz2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1))) - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][3] = dt * tz2 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * (*(*tmp1)));
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][2][4] = 0.;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][0] = dt * tz2 * (-(u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1))) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1))) + 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2)))) - dt * tz1 * (-(*(*r43)) * (*(*c34)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]);
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][1] = dt * tz2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * (*(*tmp1))));
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][2] = dt * tz2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * (*(*tmp1))));
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][3] = dt * tz2 * (2. - 0.40000000000000002) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1))) - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][3][4] = dt * tz2 * 0.40000000000000002;
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][0] = dt * tz2 * ((0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2)) - 1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][4] * (*(*tmp1)))) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1)))) - dt * tz1 * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1]))) - ((*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2]))) - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3)) * (((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]))) - (*(*c1345)) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][4]);
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][1] = dt * tz2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2))) - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1];
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][2] = dt * tz2 * (-0.40000000000000002 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2))) - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2];
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][3] = dt * tz2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][4] * (*(*tmp1))) - 0.5 * 0.40000000000000002 * ((u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][1] + u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][2] + 3. * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3]) * (*(*tmp2)))) - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2)) * u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3];
c[(i_darts_counter_temp7428)][j_darts_counter_temp7428][4][4] = dt * tz2 * (1.3999999999999999 * (u[(i_darts_counter_temp7428)][j_darts_counter_temp7428][(*(*k)) + 1][3] * (*(*tmp1)))) - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
}
(*j) = j_darts_counter_temp7428;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets7428[0].decDep();
}
TP7428::TP7428(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration, int in_lastIteration, TP7428** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),c1345_darts7428(new double*[this->numThreads]),c34_darts7428(new double*[this->numThreads]),i_darts7428(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts7428(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts7428(new int*[this->numThreads]),r43_darts7428(new double*[this->numThreads]),tmp1_darts7428(new double*[this->numThreads]),tmp2_darts7428(new double*[this->numThreads]),tmp3_darts7428(new double*[this->numThreads]), initIteration7428(in_initIteration), lastIteration7428(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets7428(new _barrierCodelets7428[1]) ,checkInCodelets7429(new _checkInCodelets7429[this->numThreads]){
/*Initialize the loop parameters*/
range7428 = abs (lastIteration7428 - initIteration7428) / 1;
rangePerCodelet7428 = range7428 / numThreads;
minIteration7428 = min<int>(lastIteration7428, initIteration7428);
remainderRange7428 = range7428 % numThreads;
/*Initialize inputs and vars.*/
this->c1345_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->c34_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->k_darts7428 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->r43_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp1_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp2_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp3_darts7428 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets7428[0] = _barrierCodelets7428(this->numThreads,this->numThreads,this, 0);
_checkInCodelets7429 * checkInCodelets7429Ptr = (this->checkInCodelets7429);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets7429);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets7429Ptr) = _checkInCodelets7429(2,1,this,codeletCounter);
#else
(*checkInCodelets7429Ptr) = _checkInCodelets7429(1,1,this,codeletCounter);
#endif
(*checkInCodelets7429Ptr).decDep();
checkInCodelets7429Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP7428::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets7429[localID].setID (codeletID);
this->checkInCodelets7429[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets7429[localID + this->baseNumThreads * i] = _checkInCodelets7429(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets7429[localID + this->baseNumThreads * i] = _checkInCodelets7429(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets7429[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets7429[localID + this->baseNumThreads * i].decDep();
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
TP7428::~TP7428(){
delete [] c1345_darts7428;
delete [] c34_darts7428;
delete [] k_darts7428;
delete [] r43_darts7428;
delete [] tmp1_darts7428;
delete [] tmp2_darts7428;
delete [] tmp3_darts7428;
delete [] barrierCodelets7428;
delete [] checkInCodelets7429;
}
/*TP9879: OMPParallelDirective*/
void TP9879::_barrierCodelets9879::fire(void)
{
TP9879* myTP = static_cast<TP9879*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP9879::_checkInCodelets9881::fire(void)
{
/*Init the vars for this region*/

/*printing node 9881: DeclStmt*/

/*printing node 9882: DeclStmt*/
this->inputsTPParent->sum0_darts9879[this->getID()] = 0.;
this->inputsTPParent->sum1_darts9879[this->getID()] = 0.;
this->inputsTPParent->sum2_darts9879[this->getID()] = 0.;
this->inputsTPParent->sum3_darts9879[this->getID()] = 0.;
this->inputsTPParent->sum4_darts9879[this->getID()] = 0.;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 9881 nextRegion: 9888 */
myTP->controlTPParent->checkInCodelets9888[this->getID()].decDep();
}
void TP9879::_checkInCodelets9888::fire(void)
{
/*Select the thread executing OMPSingleDirective 9888*/
if (! __sync_val_compare_and_swap(&(myTP->TP9888_alreadyLaunched), 0, 1))
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->sum_darts9888 = (this->inputsTPParent->sum_darts9879)/*OMP_SHARED - VAR INLINED*/;

/*printing node 9889: ForStmt*/
{
/*Loop's init*/
(this->inputsTPParent->m_darts9888) = 0;
int m_darts_counter_temp9888 = (this->inputsTPParent->m_darts9888);
for(;m_darts_counter_temp9888 < 5;m_darts_counter_temp9888++){
(*(this->inputsTPParent->sum_darts9888))[m_darts_counter_temp9888] = 0.;
}
(this->inputsTPParent->m_darts9888) = m_darts_counter_temp9888;
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp region's barrier*/
myTP->controlTPParent->barrierCodelets9888[0].decDep();
}
else
{
/*Signaling omp region's barrier*/
myTP->barrierCodelets9888[0].decDep();
}
}
void TP9879::_barrierCodelets9888::fire(void)
{
TP9879* myTP =  static_cast<TP9879*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets9897[codeletsCounter].decDep();
}
}
}
void TP9879::_checkInCodelets9897::fire(void)
{
/*region 9897 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP9897;
if(idx < myTP->TPsToUse9897){
if (!__sync_val_compare_and_swap (&(myTP->TP9897_alreadyLaunched[idx]), 0, 1)){
int range = abs ((*(this->inputsTPParent->iend_darts9879)) - (*(this->inputsTPParent->ist_darts9879))) / 1;
int rangePerCodelet = range / myTP->TPsToUse9897;
int minIteration = min<int >((*(this->inputsTPParent->iend_darts9879)), (*(this->inputsTPParent->ist_darts9879)));
int remainderRange = range % myTP->TPsToUse9897;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if((*(this->inputsTPParent->ist_darts9879)) < (*(this->inputsTPParent->iend_darts9879)))
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse9897 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse9897 - 1)
{
lastIteration = (*(this->inputsTPParent->iend_darts9879));
}
#if USEINVOKE == 1
invoke < TP9897 > (myTP, myTP->codeletsPerTP9897 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9879)), &(*(this->inputsTPParent->ist_darts9879)), &(*(this->inputsTPParent->jend_darts9879)), &(*(this->inputsTPParent->jst_darts9879)), &(*(this->inputsTPParent->nz0_darts9879)), &(myTP->TP9897Ptr[idx]));
#else
place < TP9897 > (idx, myTP, myTP->codeletsPerTP9897 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9879)), &(*(this->inputsTPParent->ist_darts9879)), &(*(this->inputsTPParent->jend_darts9879)), &(*(this->inputsTPParent->jst_darts9879)), &(*(this->inputsTPParent->nz0_darts9879)), &(myTP->TP9897Ptr[idx]));
#endif
}else{
if (myTP->TP9897Ptr[idx] != nullptr){
myTP->TP9897Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
else
{
/*Signaling next codelet region: 9897 nextRegion: 9994 */
myTP->controlTPParent->checkInCodelets9994[this->getID()].decDep();
}
}
void TP9879::_checkInCodelets9994::fire(void)
{

/*printing node 9994: CompoundAssignOperator*/
TP9992mutex.lock();
(*(this->inputsTPParent->sum_darts9879))[0] += (this->inputsTPParent->sum0_darts9879[this->getID()]);

/*printing node 9996: CompoundAssignOperator*/
(*(this->inputsTPParent->sum_darts9879))[1] += (this->inputsTPParent->sum1_darts9879[this->getID()]);

/*printing node 9998: CompoundAssignOperator*/
(*(this->inputsTPParent->sum_darts9879))[2] += (this->inputsTPParent->sum2_darts9879[this->getID()]);

/*printing node 10000: CompoundAssignOperator*/
(*(this->inputsTPParent->sum_darts9879))[3] += (this->inputsTPParent->sum3_darts9879[this->getID()]);

/*printing node 10002: CompoundAssignOperator*/
(*(this->inputsTPParent->sum_darts9879))[4] += (this->inputsTPParent->sum4_darts9879[this->getID()]);
TP9992mutex.unlock();
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 9994 nextRegion: 10004 */
myTP->controlTPParent->barrierCodelets10004[0].decDep();
}
void TP9879::_barrierCodelets10004::fire(void)
{
TP9879* myTP =  static_cast<TP9879*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets10005[codeletsCounter].decDep();
}
}
}
void TP9879::_checkInCodelets10005::fire(void)
{
/*Select the thread executing OMPSingleDirective 10005*/
if (! __sync_val_compare_and_swap(&(myTP->TP10005_alreadyLaunched), 0, 1))
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->nx0_darts10005 = (this->inputsTPParent->nx0_darts9879)/*OMP_SHARED - VAR INLINED*/;
this->inputsTPParent->ny0_darts10005 = (this->inputsTPParent->ny0_darts9879)/*OMP_SHARED - VAR INLINED*/;
this->inputsTPParent->nz0_darts10005 = (this->inputsTPParent->nz0_darts9879)/*OMP_SHARED - VAR INLINED*/;
this->inputsTPParent->sum_darts10005 = (this->inputsTPParent->sum_darts9879)/*OMP_SHARED - VAR INLINED*/;

/*printing node 10006: ForStmt*/
{
/*Loop's init*/
(this->inputsTPParent->m_darts10005) = 0;
int m_darts_counter_temp10005 = (this->inputsTPParent->m_darts10005);
for(;m_darts_counter_temp10005 < 5;m_darts_counter_temp10005++){
(*(this->inputsTPParent->sum_darts10005))[m_darts_counter_temp10005] = sqrt((*(this->inputsTPParent->sum_darts10005))[m_darts_counter_temp10005] / (((*(this->inputsTPParent->nx0_darts10005)) - 2) * ((*(this->inputsTPParent->ny0_darts10005)) - 2) * ((*(this->inputsTPParent->nz0_darts10005)) - 2)));
}
(this->inputsTPParent->m_darts10005) = m_darts_counter_temp10005;
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp region's barrier*/
myTP->controlTPParent->barrierCodelets10005[0].decDep();
}
else
{
/*Signaling omp region's barrier*/
myTP->barrierCodelets10005[0].decDep();
}
}
void TP9879::_barrierCodelets10005::fire(void)
{
TP9879* myTP =  static_cast<TP9879*>(myTP_);
myTP->TPParent->barrierCodelets9879[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets9879[0]));
}
TP9879::TP9879(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int *in_iend, int *in_ist, int *in_jend, int *in_jst, int *in_nx0, int *in_ny0, int *in_nz0, double * *in_sum):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),iend_darts9879(in_iend)/*OMP_SHARED - INPUT*/,ist_darts9879(in_ist)/*OMP_SHARED - INPUT*/,jend_darts9879(in_jend)/*OMP_SHARED - INPUT*/,jst_darts9879(in_jst)/*OMP_SHARED - INPUT*/,nx0_darts9879(in_nx0)/*OMP_SHARED - INPUT*/,ny0_darts9879(in_ny0)/*OMP_SHARED - INPUT*/,nz0_darts9879(in_nz0)/*OMP_SHARED - INPUT*/,sum_darts9879(in_sum)/*OMP_SHARED - INPUT*/,i_darts9879(new int[this->numThreads])/*VARIABLE*/,j_darts9879(new int[this->numThreads])/*VARIABLE*/,k_darts9879(new int[this->numThreads])/*VARIABLE*/,m_darts9879(new int[this->numThreads])/*VARIABLE*/,sum0_darts9879(new double[this->numThreads])/*VARIABLE*/,sum1_darts9879(new double[this->numThreads])/*VARIABLE*/,sum2_darts9879(new double[this->numThreads])/*VARIABLE*/,sum3_darts9879(new double[this->numThreads])/*VARIABLE*/,sum4_darts9879(new double[this->numThreads])/*VARIABLE*/, TP9888_alreadyLaunched(0), TP9897Ptr(new TP9897 *[NUMTPS9897]), TP9897_alreadyLaunched(new size_t [NUMTPS9897]), numTPsSet9897(0), numTPsReady9897(0), TPsToUse9897(NUMTPS9897), codeletsPerTP9897(this->numThreads/NUMTPS9897), totalCodelets9897(this->TPsToUse9897*this->codeletsPerTP9897), TP10005_alreadyLaunched(0) ,barrierCodelets9879(new _barrierCodelets9879[1]) ,checkInCodelets9881(new _checkInCodelets9881[this->numThreads]) ,checkInCodelets9888(new _checkInCodelets9888[this->numThreads]) ,barrierCodelets9888(new _barrierCodelets9888[1]) ,checkInCodelets9897(new _checkInCodelets9897[this->numThreads]) ,checkInCodelets9994(new _checkInCodelets9994[this->numThreads]) ,barrierCodelets10004(new _barrierCodelets10004[1]) ,checkInCodelets10005(new _checkInCodelets10005[this->numThreads]) ,barrierCodelets10005(new _barrierCodelets10005[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets9879[0] = _barrierCodelets9879(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets10005[0] = _barrierCodelets10005(this->numThreads,this->numThreads,this, 0);
barrierCodelets10004[0] = _barrierCodelets10004(this->numThreads,this->numThreads,this, 0);
barrierCodelets9888[0] = _barrierCodelets9888(this->numThreads,this->numThreads,this, 0);
_checkInCodelets10005 * checkInCodelets10005Ptr = (this->checkInCodelets10005);
_checkInCodelets9994 * checkInCodelets9994Ptr = (this->checkInCodelets9994);
_checkInCodelets9897 * checkInCodelets9897Ptr = (this->checkInCodelets9897);
for(int i=0; i<NUMTPS9897; i++)
{
TP9897Ptr[i] = nullptr;
TP9897_alreadyLaunched[i] = 0;
}
_checkInCodelets9888 * checkInCodelets9888Ptr = (this->checkInCodelets9888);
_checkInCodelets9881 * checkInCodelets9881Ptr = (this->checkInCodelets9881);
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets10005Ptr) = _checkInCodelets10005(1,1,this,codeletCounter);
checkInCodelets10005Ptr++;
(*checkInCodelets9994Ptr) = _checkInCodelets9994(1,1,this,codeletCounter);
checkInCodelets9994Ptr++;
(*checkInCodelets9897Ptr) = _checkInCodelets9897(1,1,this,codeletCounter);
checkInCodelets9897Ptr++;
(*checkInCodelets9888Ptr) = _checkInCodelets9888(1,1,this,codeletCounter);
checkInCodelets9888Ptr++;
(*checkInCodelets9881Ptr) = _checkInCodelets9881(1,1,this,codeletCounter);
(*checkInCodelets9881Ptr).decDep();
checkInCodelets9881Ptr++;
}
}
TP9879::~TP9879(){
delete []i_darts9879;
delete []j_darts9879;
delete []k_darts9879;
delete []m_darts9879;
delete []sum0_darts9879;
delete []sum1_darts9879;
delete []sum2_darts9879;
delete []sum3_darts9879;
delete []sum4_darts9879;
delete [] barrierCodelets9879;
delete [] barrierCodelets10005;
delete [] checkInCodelets10005;
delete [] barrierCodelets10004;
delete [] checkInCodelets9994;
delete [] checkInCodelets9897;
delete [] barrierCodelets9888;
delete [] checkInCodelets9888;
delete [] checkInCodelets9881;
}
/*TP9897: OMPForDirective*/
bool TP9897::requestNewRangeIterations9897(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet9897 * codeletID;
int tempEndRange   = rangePerCodelet9897 * (codeletID + 1);
if (remainderRange9897 != 0)
{
if (codeletID < (uint32_t)remainderRange9897)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange9897;
tempEndRange += remainderRange9897;
}
}
tempStartRange = tempStartRange*1 + minIteration9897;
tempEndRange = tempEndRange*1 + minIteration9897;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration9897 < lastIteration9897)
{
(this->inputsTPParent->i_darts9897[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts9897[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration9897;
}
}
return isThereNewIteration;
}
void TP9897::_checkInCodelets9898::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->sum0_darts9897[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->sum0_darts9879[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->sum1_darts9897[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->sum1_darts9879[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->sum2_darts9897[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->sum2_darts9879[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->sum3_darts9897[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->sum3_darts9879[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->sum4_darts9897[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->sum4_darts9879[this->getID()]);

/*printing node 9898: ForStmt*/
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
int* i = &(this->inputsTPParent->i_darts9897[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts9897[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* jend = (this->inputsTPParent->jend_darts9897);
(void)jend/*OMP_SHARED*/;
int* jst = (this->inputsTPParent->jst_darts9897);
(void)jst/*OMP_SHARED*/;
int* k = &(this->inputsTPParent->k_darts9897[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* nz0 = (this->inputsTPParent->nz0_darts9897);
(void)nz0/*OMP_SHARED*/;
double** sum0 = &(this->inputsTPParent->sum0_darts9897[this->getLocalID()]);
(void)sum0/*OMP_SHARED_PRIVATE*/;
double** sum1 = &(this->inputsTPParent->sum1_darts9897[this->getLocalID()]);
(void)sum1/*OMP_SHARED_PRIVATE*/;
double** sum2 = &(this->inputsTPParent->sum2_darts9897[this->getLocalID()]);
(void)sum2/*OMP_SHARED_PRIVATE*/;
double** sum3 = &(this->inputsTPParent->sum3_darts9897[this->getLocalID()]);
(void)sum3/*OMP_SHARED_PRIVATE*/;
double** sum4 = &(this->inputsTPParent->sum4_darts9897[this->getLocalID()]);
(void)sum4/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations9897((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Find and signal the next codelet*/
myTP->controlTPParent->TPParent->checkInCodelets9994[this->getID()].decDep();
return;
}
for (int i_darts_counter_temp9897 = (*i);i_darts_counter_temp9897<=endRange && i_darts_counter_temp9897<=this->inputsTPParent->lastIteration9897;i_darts_counter_temp9897++)
{
{
{
/*Loop's init*/
(*j) = (*(jst));
int j_darts_counter_temp9897 = (*j);
for(;j_darts_counter_temp9897 <= (*(jend));j_darts_counter_temp9897++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp9897 = (*k);
for(;k_darts_counter_temp9897 <= (*(nz0)) - 2;k_darts_counter_temp9897++){
(*(*sum0)) = (*(*sum0)) + rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][0] * rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][0];
(*(*sum1)) = (*(*sum1)) + rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][1] * rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][1];
(*(*sum2)) = (*(*sum2)) + rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][2] * rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][2];
(*(*sum3)) = (*(*sum3)) + rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][3] * rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][3];
(*(*sum4)) = (*(*sum4)) + rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][4] * rsd[(i_darts_counter_temp9897)][j_darts_counter_temp9897][k_darts_counter_temp9897][4];
}
(*k) = k_darts_counter_temp9897;
}
}
(*j) = j_darts_counter_temp9897;
}
}
}
/*If this omp for has no barrier, 
check if all the codelets 
replicated from the same 
global ID has finished and 
signal the next codelet. 
Otherwise, return.*/
uint32_t completedMultCodelet = __sync_fetch_and_add(&(myTP->signalNextReady[this->getLocalID() % myTP->baseNumThreads]), 1);
if(completedMultCodelet < (uint32_t)(DARTS_CODELETS_MULT - 1))
return;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/
myTP->controlTPParent->TPParent->checkInCodelets9994[this->getID()].decDep();
}
TP9897::TP9897(int in_numThreads, int in_mainCodeletID, TP9879* in_TPParent, int in_initIteration, int in_lastIteration, int *in_iend, int *in_ist, int *in_jend, int *in_jst, int *in_nz0, TP9897** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts9897(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iend_darts9897(in_iend)/*OMP_SHARED - INPUT*/,ist_darts9897(in_ist)/*OMP_SHARED - INPUT*/,j_darts9897(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jend_darts9897(in_jend)/*OMP_SHARED - INPUT*/,jst_darts9897(in_jst)/*OMP_SHARED - INPUT*/,k_darts9897(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,nz0_darts9897(in_nz0)/*OMP_SHARED - INPUT*/,sum0_darts9897(new double*[this->numThreads]),sum1_darts9897(new double*[this->numThreads]),sum2_darts9897(new double*[this->numThreads]),sum3_darts9897(new double*[this->numThreads]),sum4_darts9897(new double*[this->numThreads]), initIteration9897(in_initIteration), lastIteration9897(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT), signalNextReady(new int[baseNumThreads]) ,checkInCodelets9898(new _checkInCodelets9898[this->numThreads]){
/*Initialize the loop parameters*/
range9897 = abs (lastIteration9897 - initIteration9897) / 1;
rangePerCodelet9897 = range9897 / numThreads;
minIteration9897 = min<int>(lastIteration9897, initIteration9897);
remainderRange9897 = range9897 % numThreads;
/*Initialize inputs and vars.*/
this->sum0_darts9897 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->sum1_darts9897 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->sum2_darts9897 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->sum3_darts9897 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->sum4_darts9897 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
_checkInCodelets9898 * checkInCodelets9898Ptr = (this->checkInCodelets9898);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets9898);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets9898Ptr) = _checkInCodelets9898(2,1,this,codeletCounter);
#else
(*checkInCodelets9898Ptr) = _checkInCodelets9898(1,1,this,codeletCounter);
#endif
(*checkInCodelets9898Ptr).decDep();
checkInCodelets9898Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP9897::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets9898[localID].setID (codeletID);
this->checkInCodelets9898[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets9898[localID + this->baseNumThreads * i] = _checkInCodelets9898(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets9898[localID + this->baseNumThreads * i] = _checkInCodelets9898(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets9898[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets9898[localID + this->baseNumThreads * i].decDep();
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
TP9897::~TP9897(){
delete [] sum0_darts9897;
delete [] sum1_darts9897;
delete [] sum2_darts9897;
delete [] sum3_darts9897;
delete [] sum4_darts9897;
delete [] checkInCodelets9898;
}
/*TP10796: OMPParallelDirective*/
void TP10796::_barrierCodelets10796::fire(void)
{
TP10796* myTP = static_cast<TP10796*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP10796::_checkInCodelets10811::fire(void)
{
/*region 10811 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP10811;
if(idx < myTP->TPsToUse10811){
if (!__sync_val_compare_and_swap (&(myTP->TP10811_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 1 - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse10811;
int minIteration = min<int >(nx - 1, 0);
int remainderRange = range % myTP->TPsToUse10811;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx - 1)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse10811 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse10811 - 1)
{
lastIteration = nx - 1;
}
#if USEINVOKE == 1
invoke < TP10811 > (myTP, myTP->codeletsPerTP10811 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP10811Ptr[idx]));
#else
place < TP10811 > (idx, myTP, myTP->codeletsPerTP10811 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP10811Ptr[idx]));
#endif
}else{
if (myTP->TP10811Ptr[idx] != nullptr){
myTP->TP10811Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets10811::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets10868[codeletsCounter].decDep();
}
}
}
void TP10796::_checkInCodelets10868::fire(void)
{

/*printing node 10868: BinaryOperator*/
(this->inputsTPParent->L1_darts10796[this->getID()]) = 0;

/*printing node 10869: BinaryOperator*/
(this->inputsTPParent->L2_darts10796[this->getID()]) = nx - 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 10868 nextRegion: 10871 */
myTP->controlTPParent->checkInCodelets10871[this->getID()].decDep();
}
void TP10796::_checkInCodelets10871::fire(void)
{
/*region 10871 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP10871;
if(idx < myTP->TPsToUse10871){
if (!__sync_val_compare_and_swap (&(myTP->TP10871_alreadyLaunched[idx]), 0, 1)){
int range = abs ((this->inputsTPParent->L2_darts10796[this->getID()]) - (this->inputsTPParent->L1_darts10796[this->getID()])) / 1;
int rangePerCodelet = range / myTP->TPsToUse10871;
int minIteration = min<int >((this->inputsTPParent->L2_darts10796[this->getID()]), (this->inputsTPParent->L1_darts10796[this->getID()]));
int remainderRange = range % myTP->TPsToUse10871;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if((this->inputsTPParent->L1_darts10796[this->getID()]) < (this->inputsTPParent->L2_darts10796[this->getID()]))
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse10871 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse10871 - 1)
{
lastIteration = (this->inputsTPParent->L2_darts10796[this->getID()]);
}
#if USEINVOKE == 1
invoke < TP10871 > (myTP, myTP->codeletsPerTP10871 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP10871Ptr[idx]));
#else
place < TP10871 > (idx, myTP, myTP->codeletsPerTP10871 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP10871Ptr[idx]));
#endif
}else{
if (myTP->TP10871Ptr[idx] != nullptr){
myTP->TP10871Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets10871::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets11020[codeletsCounter].decDep();
}
}
}
void TP10796::_checkInCodelets11020::fire(void)
{
/*region 11020 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP11020;
if(idx < myTP->TPsToUse11020){
if (!__sync_val_compare_and_swap (&(myTP->TP11020_alreadyLaunched[idx]), 0, 1)){
int range = abs (jend - jst) / 1;
int rangePerCodelet = range / myTP->TPsToUse11020;
int minIteration = min<int >(jend, jst);
int remainderRange = range % myTP->TPsToUse11020;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(jst < jend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse11020 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse11020 - 1)
{
lastIteration = jend;
}
#if USEINVOKE == 1
invoke < TP11020 > (myTP, myTP->codeletsPerTP11020 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11020Ptr[idx]));
#else
place < TP11020 > (idx, myTP, myTP->codeletsPerTP11020 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11020Ptr[idx]));
#endif
}else{
if (myTP->TP11020Ptr[idx] != nullptr){
myTP->TP11020Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets11020::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets11657[codeletsCounter].decDep();
}
}
}
void TP10796::_checkInCodelets11657::fire(void)
{

/*printing node 11657: BinaryOperator*/
(this->inputsTPParent->L1_darts10796[this->getID()]) = 0;

/*printing node 11658: BinaryOperator*/
(this->inputsTPParent->L2_darts10796[this->getID()]) = ny - 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 11657 nextRegion: 11660 */
myTP->controlTPParent->checkInCodelets11660[this->getID()].decDep();
}
void TP10796::_checkInCodelets11660::fire(void)
{
/*region 11660 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP11660;
if(idx < myTP->TPsToUse11660){
if (!__sync_val_compare_and_swap (&(myTP->TP11660_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse11660;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse11660;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse11660 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse11660 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP11660 > (myTP, myTP->codeletsPerTP11660 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11660Ptr[idx]));
#else
place < TP11660 > (idx, myTP, myTP->codeletsPerTP11660 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11660Ptr[idx]));
#endif
}else{
if (myTP->TP11660Ptr[idx] != nullptr){
myTP->TP11660Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets11660::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets11809[codeletsCounter].decDep();
}
}
}
void TP10796::_checkInCodelets11809::fire(void)
{
/*region 11809 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP11809;
if(idx < myTP->TPsToUse11809){
if (!__sync_val_compare_and_swap (&(myTP->TP11809_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse11809;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse11809;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse11809 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse11809 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP11809 > (myTP, myTP->codeletsPerTP11809 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11809Ptr[idx]));
#else
place < TP11809 > (idx, myTP, myTP->codeletsPerTP11809 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP11809Ptr[idx]));
#endif
}else{
if (myTP->TP11809Ptr[idx] != nullptr){
myTP->TP11809Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets11809::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets12446[codeletsCounter].decDep();
}
}
}
void TP10796::_checkInCodelets12446::fire(void)
{
/*region 12446 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP12446;
if(idx < myTP->TPsToUse12446){
if (!__sync_val_compare_and_swap (&(myTP->TP12446_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse12446;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse12446;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse12446 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse12446 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP12446 > (myTP, myTP->codeletsPerTP12446 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP12446Ptr[idx]));
#else
place < TP12446 > (idx, myTP, myTP->codeletsPerTP12446 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP12446Ptr[idx]));
#endif
}else{
if (myTP->TP12446Ptr[idx] != nullptr){
myTP->TP12446Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP10796::_barrierCodelets12446::fire(void)
{
TP10796* myTP =  static_cast<TP10796*>(myTP_);
myTP->TPParent->barrierCodelets10796[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets10796[0]));
}
TP10796::TP10796(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),L1_darts10796(new int[this->numThreads])/*VARIABLE*/,L2_darts10796(new int[this->numThreads])/*VARIABLE*/,i_darts10796(new int[this->numThreads])/*VARIABLE*/,iend1_darts10796(new int[this->numThreads])/*VARIABLE*/,ist1_darts10796(new int[this->numThreads])/*VARIABLE*/,j_darts10796(new int[this->numThreads])/*VARIABLE*/,jend1_darts10796(new int[this->numThreads])/*VARIABLE*/,jst1_darts10796(new int[this->numThreads])/*VARIABLE*/,k_darts10796(new int[this->numThreads])/*VARIABLE*/,m_darts10796(new int[this->numThreads])/*VARIABLE*/,q_darts10796(new double[this->numThreads])/*VARIABLE*/,tmp_darts10796(new double[this->numThreads])/*VARIABLE*/,u21_darts10796(new double[this->numThreads])/*VARIABLE*/,u21i_darts10796(new double[this->numThreads])/*VARIABLE*/,u21im1_darts10796(new double[this->numThreads])/*VARIABLE*/,u21j_darts10796(new double[this->numThreads])/*VARIABLE*/,u21jm1_darts10796(new double[this->numThreads])/*VARIABLE*/,u21k_darts10796(new double[this->numThreads])/*VARIABLE*/,u21km1_darts10796(new double[this->numThreads])/*VARIABLE*/,u31_darts10796(new double[this->numThreads])/*VARIABLE*/,u31i_darts10796(new double[this->numThreads])/*VARIABLE*/,u31im1_darts10796(new double[this->numThreads])/*VARIABLE*/,u31j_darts10796(new double[this->numThreads])/*VARIABLE*/,u31jm1_darts10796(new double[this->numThreads])/*VARIABLE*/,u31k_darts10796(new double[this->numThreads])/*VARIABLE*/,u31km1_darts10796(new double[this->numThreads])/*VARIABLE*/,u41_darts10796(new double[this->numThreads])/*VARIABLE*/,u41i_darts10796(new double[this->numThreads])/*VARIABLE*/,u41im1_darts10796(new double[this->numThreads])/*VARIABLE*/,u41j_darts10796(new double[this->numThreads])/*VARIABLE*/,u41jm1_darts10796(new double[this->numThreads])/*VARIABLE*/,u41k_darts10796(new double[this->numThreads])/*VARIABLE*/,u41km1_darts10796(new double[this->numThreads])/*VARIABLE*/,u51i_darts10796(new double[this->numThreads])/*VARIABLE*/,u51im1_darts10796(new double[this->numThreads])/*VARIABLE*/,u51j_darts10796(new double[this->numThreads])/*VARIABLE*/,u51jm1_darts10796(new double[this->numThreads])/*VARIABLE*/,u51k_darts10796(new double[this->numThreads])/*VARIABLE*/,u51km1_darts10796(new double[this->numThreads])/*VARIABLE*/, TP10811Ptr(new TP10811 *[NUMTPS10811]), TP10811_alreadyLaunched(new size_t [NUMTPS10811]), numTPsSet10811(0), numTPsReady10811(0), TPsToUse10811(NUMTPS10811), codeletsPerTP10811(this->numThreads/NUMTPS10811), totalCodelets10811(this->TPsToUse10811*this->codeletsPerTP10811), TP10871Ptr(new TP10871 *[NUMTPS10871]), TP10871_alreadyLaunched(new size_t [NUMTPS10871]), numTPsSet10871(0), numTPsReady10871(0), TPsToUse10871(NUMTPS10871), codeletsPerTP10871(this->numThreads/NUMTPS10871), totalCodelets10871(this->TPsToUse10871*this->codeletsPerTP10871), TP11020Ptr(new TP11020 *[NUMTPS11020]), TP11020_alreadyLaunched(new size_t [NUMTPS11020]), numTPsSet11020(0), numTPsReady11020(0), TPsToUse11020(NUMTPS11020), codeletsPerTP11020(this->numThreads/NUMTPS11020), totalCodelets11020(this->TPsToUse11020*this->codeletsPerTP11020), TP11660Ptr(new TP11660 *[NUMTPS11660]), TP11660_alreadyLaunched(new size_t [NUMTPS11660]), numTPsSet11660(0), numTPsReady11660(0), TPsToUse11660(NUMTPS11660), codeletsPerTP11660(this->numThreads/NUMTPS11660), totalCodelets11660(this->TPsToUse11660*this->codeletsPerTP11660), TP11809Ptr(new TP11809 *[NUMTPS11809]), TP11809_alreadyLaunched(new size_t [NUMTPS11809]), numTPsSet11809(0), numTPsReady11809(0), TPsToUse11809(NUMTPS11809), codeletsPerTP11809(this->numThreads/NUMTPS11809), totalCodelets11809(this->TPsToUse11809*this->codeletsPerTP11809), TP12446Ptr(new TP12446 *[NUMTPS12446]), TP12446_alreadyLaunched(new size_t [NUMTPS12446]), numTPsSet12446(0), numTPsReady12446(0), TPsToUse12446(NUMTPS12446), codeletsPerTP12446(this->numThreads/NUMTPS12446), totalCodelets12446(this->TPsToUse12446*this->codeletsPerTP12446) ,barrierCodelets10796(new _barrierCodelets10796[1]) ,checkInCodelets10811(new _checkInCodelets10811[this->numThreads]) ,barrierCodelets10811(new _barrierCodelets10811[1]) ,checkInCodelets10868(new _checkInCodelets10868[this->numThreads]) ,checkInCodelets10871(new _checkInCodelets10871[this->numThreads]) ,barrierCodelets10871(new _barrierCodelets10871[1]) ,checkInCodelets11020(new _checkInCodelets11020[this->numThreads]) ,barrierCodelets11020(new _barrierCodelets11020[1]) ,checkInCodelets11657(new _checkInCodelets11657[this->numThreads]) ,checkInCodelets11660(new _checkInCodelets11660[this->numThreads]) ,barrierCodelets11660(new _barrierCodelets11660[1]) ,checkInCodelets11809(new _checkInCodelets11809[this->numThreads]) ,barrierCodelets11809(new _barrierCodelets11809[1]) ,checkInCodelets12446(new _checkInCodelets12446[this->numThreads]) ,barrierCodelets12446(new _barrierCodelets12446[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets10796[0] = _barrierCodelets10796(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets12446[0] = _barrierCodelets12446(NUMTPS12446,NUMTPS12446,this, 0);
barrierCodelets11809[0] = _barrierCodelets11809(NUMTPS11809,NUMTPS11809,this, 0);
barrierCodelets11660[0] = _barrierCodelets11660(NUMTPS11660,NUMTPS11660,this, 0);
barrierCodelets11020[0] = _barrierCodelets11020(NUMTPS11020,NUMTPS11020,this, 0);
barrierCodelets10871[0] = _barrierCodelets10871(NUMTPS10871,NUMTPS10871,this, 0);
barrierCodelets10811[0] = _barrierCodelets10811(NUMTPS10811,NUMTPS10811,this, 0);
_checkInCodelets12446 * checkInCodelets12446Ptr = (this->checkInCodelets12446);
for(int i=0; i<NUMTPS12446; i++)
{
TP12446Ptr[i] = nullptr;
TP12446_alreadyLaunched[i] = 0;
}
_checkInCodelets11809 * checkInCodelets11809Ptr = (this->checkInCodelets11809);
for(int i=0; i<NUMTPS11809; i++)
{
TP11809Ptr[i] = nullptr;
TP11809_alreadyLaunched[i] = 0;
}
_checkInCodelets11660 * checkInCodelets11660Ptr = (this->checkInCodelets11660);
for(int i=0; i<NUMTPS11660; i++)
{
TP11660Ptr[i] = nullptr;
TP11660_alreadyLaunched[i] = 0;
}
_checkInCodelets11657 * checkInCodelets11657Ptr = (this->checkInCodelets11657);
_checkInCodelets11020 * checkInCodelets11020Ptr = (this->checkInCodelets11020);
for(int i=0; i<NUMTPS11020; i++)
{
TP11020Ptr[i] = nullptr;
TP11020_alreadyLaunched[i] = 0;
}
_checkInCodelets10871 * checkInCodelets10871Ptr = (this->checkInCodelets10871);
for(int i=0; i<NUMTPS10871; i++)
{
TP10871Ptr[i] = nullptr;
TP10871_alreadyLaunched[i] = 0;
}
_checkInCodelets10868 * checkInCodelets10868Ptr = (this->checkInCodelets10868);
_checkInCodelets10811 * checkInCodelets10811Ptr = (this->checkInCodelets10811);
for(int i=0; i<NUMTPS10811; i++)
{
TP10811Ptr[i] = nullptr;
TP10811_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets12446Ptr) = _checkInCodelets12446(1,1,this,codeletCounter);
checkInCodelets12446Ptr++;
(*checkInCodelets11809Ptr) = _checkInCodelets11809(1,1,this,codeletCounter);
checkInCodelets11809Ptr++;
(*checkInCodelets11660Ptr) = _checkInCodelets11660(1,1,this,codeletCounter);
checkInCodelets11660Ptr++;
(*checkInCodelets11657Ptr) = _checkInCodelets11657(1,1,this,codeletCounter);
checkInCodelets11657Ptr++;
(*checkInCodelets11020Ptr) = _checkInCodelets11020(1,1,this,codeletCounter);
checkInCodelets11020Ptr++;
(*checkInCodelets10871Ptr) = _checkInCodelets10871(1,1,this,codeletCounter);
checkInCodelets10871Ptr++;
(*checkInCodelets10868Ptr) = _checkInCodelets10868(1,1,this,codeletCounter);
checkInCodelets10868Ptr++;
(*checkInCodelets10811Ptr) = _checkInCodelets10811(1,1,this,codeletCounter);
(*checkInCodelets10811Ptr).decDep();
checkInCodelets10811Ptr++;
}
}
TP10796::~TP10796(){
delete []L1_darts10796;
delete []L2_darts10796;
delete []i_darts10796;
delete []iend1_darts10796;
delete []ist1_darts10796;
delete []j_darts10796;
delete []jend1_darts10796;
delete []jst1_darts10796;
delete []k_darts10796;
delete []m_darts10796;
delete []q_darts10796;
delete []tmp_darts10796;
delete []u21_darts10796;
delete []u21i_darts10796;
delete []u21im1_darts10796;
delete []u21j_darts10796;
delete []u21jm1_darts10796;
delete []u21k_darts10796;
delete []u21km1_darts10796;
delete []u31_darts10796;
delete []u31i_darts10796;
delete []u31im1_darts10796;
delete []u31j_darts10796;
delete []u31jm1_darts10796;
delete []u31k_darts10796;
delete []u31km1_darts10796;
delete []u41_darts10796;
delete []u41i_darts10796;
delete []u41im1_darts10796;
delete []u41j_darts10796;
delete []u41jm1_darts10796;
delete []u41k_darts10796;
delete []u41km1_darts10796;
delete []u51i_darts10796;
delete []u51im1_darts10796;
delete []u51j_darts10796;
delete []u51jm1_darts10796;
delete []u51k_darts10796;
delete []u51km1_darts10796;
delete [] barrierCodelets10796;
delete [] barrierCodelets12446;
delete [] checkInCodelets12446;
delete [] barrierCodelets11809;
delete [] checkInCodelets11809;
delete [] barrierCodelets11660;
delete [] checkInCodelets11660;
delete [] checkInCodelets11657;
delete [] barrierCodelets11020;
delete [] checkInCodelets11020;
delete [] barrierCodelets10871;
delete [] checkInCodelets10871;
delete [] checkInCodelets10868;
delete [] barrierCodelets10811;
delete [] checkInCodelets10811;
}
/*TP10811: OMPForDirective*/
void TP10811::_barrierCodelets10811::fire(void)
{
TP10811* myTP = static_cast<TP10811*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets10811[0].decDep ();
}
bool TP10811::requestNewRangeIterations10811(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet10811 * codeletID;
int tempEndRange   = rangePerCodelet10811 * (codeletID + 1);
if (remainderRange10811 != 0)
{
if (codeletID < (uint32_t)remainderRange10811)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange10811;
tempEndRange += remainderRange10811;
}
}
tempStartRange = tempStartRange*1 + minIteration10811;
tempEndRange = tempEndRange*1 + minIteration10811;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration10811 < lastIteration10811)
{
(this->inputsTPParent->i_darts10811[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts10811[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration10811;
}
}
return isThereNewIteration;
}
void TP10811::_checkInCodelets10812::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 10812: ForStmt*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: m*/
int* i = &(this->inputsTPParent->i_darts10811[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts10811[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts10811[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts10811[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10811((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets10811[0].decDep();
return;
}
for (int i_darts_counter_temp10811 = (*i);i_darts_counter_temp10811<=endRange && i_darts_counter_temp10811<=this->inputsTPParent->lastIteration10811;i_darts_counter_temp10811++)
{
{
{
/*Loop's init*/
(*j) = 0;
int j_darts_counter_temp10811 = (*j);
for(;j_darts_counter_temp10811 <= ny - 1;j_darts_counter_temp10811++){
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp10811 = (*k);
for(;k_darts_counter_temp10811 <= nz - 1;k_darts_counter_temp10811++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp10811 = (*m);
for(;m_darts_counter_temp10811 < 5;m_darts_counter_temp10811++){
rsd[(i_darts_counter_temp10811)][j_darts_counter_temp10811][k_darts_counter_temp10811][m_darts_counter_temp10811] = -frct[(i_darts_counter_temp10811)][j_darts_counter_temp10811][k_darts_counter_temp10811][m_darts_counter_temp10811];
}
(*m) = m_darts_counter_temp10811;
}
}
(*k) = k_darts_counter_temp10811;
}
}
(*j) = j_darts_counter_temp10811;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets10811[0].decDep();
}
TP10811::TP10811(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP10811** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts10811(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts10811(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts10811(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts10811(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration10811(in_initIteration), lastIteration10811(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets10811(new _barrierCodelets10811[1]) ,checkInCodelets10812(new _checkInCodelets10812[this->numThreads]){
/*Initialize the loop parameters*/
range10811 = abs (lastIteration10811 - initIteration10811) / 1;
rangePerCodelet10811 = range10811 / numThreads;
minIteration10811 = min<int>(lastIteration10811, initIteration10811);
remainderRange10811 = range10811 % numThreads;
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets10811[0] = _barrierCodelets10811(this->numThreads,this->numThreads,this, 0);
_checkInCodelets10812 * checkInCodelets10812Ptr = (this->checkInCodelets10812);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets10812);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets10812Ptr) = _checkInCodelets10812(2,1,this,codeletCounter);
#else
(*checkInCodelets10812Ptr) = _checkInCodelets10812(1,1,this,codeletCounter);
#endif
(*checkInCodelets10812Ptr).decDep();
checkInCodelets10812Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP10811::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets10812[localID].setID (codeletID);
this->checkInCodelets10812[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets10812[localID + this->baseNumThreads * i] = _checkInCodelets10812(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets10812[localID + this->baseNumThreads * i] = _checkInCodelets10812(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets10812[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets10812[localID + this->baseNumThreads * i].decDep();
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
TP10811::~TP10811(){
delete [] barrierCodelets10811;
delete [] checkInCodelets10812;
}
/*TP10871: OMPForDirective*/
void TP10871::_barrierCodelets10871::fire(void)
{
TP10871* myTP = static_cast<TP10871*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets10871[0].decDep ();
}
bool TP10871::requestNewRangeIterations10871(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet10871 * codeletID;
int tempEndRange   = rangePerCodelet10871 * (codeletID + 1);
if (remainderRange10871 != 0)
{
if (codeletID < (uint32_t)remainderRange10871)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange10871;
tempEndRange += remainderRange10871;
}
}
tempStartRange = tempStartRange*1 + minIteration10871;
tempEndRange = tempEndRange*1 + minIteration10871;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration10871 < lastIteration10871)
{
(this->inputsTPParent->i_darts10871[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts10871[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration10871;
}
}
return isThereNewIteration;
}
void TP10871::_checkInCodelets10872::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L1_darts10871[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts10871[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts10871[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21_darts10871[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21_darts10796[this->getID()]);

/*printing node 10872: ForStmt*/
/*var: L1*/
/*var: L2*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: q*/
/*var: u21*/
int* i = &(this->inputsTPParent->i_darts10871[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts10871[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts10871[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts10871[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** u21 = &(this->inputsTPParent->u21_darts10871[this->getLocalID()]);
(void)u21/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10871((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets10871[0].decDep();
return;
}
for (int i_darts_counter_temp10871 = (*i);i_darts_counter_temp10871<=endRange && i_darts_counter_temp10871<=this->inputsTPParent->lastIteration10871;i_darts_counter_temp10871++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp10871 = (*j);
for(;j_darts_counter_temp10871 <= jend;j_darts_counter_temp10871++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp10871 = (*k);
for(;k_darts_counter_temp10871 <= nz - 2;k_darts_counter_temp10871++){
flux[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][0] = u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1];
(*(*u21)) = u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1] / u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][0];
(*(*q)) = 0.5 * (u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1] * u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1] + u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][2] * u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][2] + u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][3] * u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][3]) / u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][0];
flux[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1] = u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][1] * (*(*u21)) + 0.40000000000000002 * (u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][4] - (*(*q)));
flux[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][2] = u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][2] * (*(*u21));
flux[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][3] = u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][3] * (*(*u21));
flux[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][4] = (1.3999999999999999 * u[(i_darts_counter_temp10871)][j_darts_counter_temp10871][k_darts_counter_temp10871][4] - 0.40000000000000002 * (*(*q))) * (*(*u21));
}
(*k) = k_darts_counter_temp10871;
}
}
(*j) = j_darts_counter_temp10871;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets10871[0].decDep();
}
TP10871::TP10871(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP10871** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L1_darts10871(new int*[this->numThreads]),L2_darts10871(new int*[this->numThreads]),i_darts10871(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts10871(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts10871(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts10871(new double*[this->numThreads]),u21_darts10871(new double*[this->numThreads]), initIteration10871(in_initIteration), lastIteration10871(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets10871(new _barrierCodelets10871[1]) ,checkInCodelets10872(new _checkInCodelets10872[this->numThreads]){
/*Initialize the loop parameters*/
range10871 = abs (lastIteration10871 - initIteration10871) / 1;
rangePerCodelet10871 = range10871 / numThreads;
minIteration10871 = min<int>(lastIteration10871, initIteration10871);
remainderRange10871 = range10871 % numThreads;
/*Initialize inputs and vars.*/
this->L1_darts10871 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->L2_darts10871 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->q_darts10871 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21_darts10871 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets10871[0] = _barrierCodelets10871(this->numThreads,this->numThreads,this, 0);
_checkInCodelets10872 * checkInCodelets10872Ptr = (this->checkInCodelets10872);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets10872);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets10872Ptr) = _checkInCodelets10872(2,1,this,codeletCounter);
#else
(*checkInCodelets10872Ptr) = _checkInCodelets10872(1,1,this,codeletCounter);
#endif
(*checkInCodelets10872Ptr).decDep();
checkInCodelets10872Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP10871::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets10872[localID].setID (codeletID);
this->checkInCodelets10872[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets10872[localID + this->baseNumThreads * i] = _checkInCodelets10872(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets10872[localID + this->baseNumThreads * i] = _checkInCodelets10872(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets10872[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets10872[localID + this->baseNumThreads * i].decDep();
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
TP10871::~TP10871(){
delete [] L1_darts10871;
delete [] L2_darts10871;
delete [] q_darts10871;
delete [] u21_darts10871;
delete [] barrierCodelets10871;
delete [] checkInCodelets10872;
}
/*TP11020: OMPForDirective*/
void TP11020::_barrierCodelets11020::fire(void)
{
TP11020* myTP = static_cast<TP11020*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets11020[0].decDep ();
}
bool TP11020::requestNewRangeIterations11020(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet11020 * codeletID;
int tempEndRange   = rangePerCodelet11020 * (codeletID + 1);
if (remainderRange11020 != 0)
{
if (codeletID < (uint32_t)remainderRange11020)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange11020;
tempEndRange += remainderRange11020;
}
}
tempStartRange = tempStartRange*1 + minIteration11020;
tempEndRange = tempEndRange*1 + minIteration11020;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration11020 < lastIteration11020)
{
(this->inputsTPParent->j_darts11020[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->j_darts11020[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration11020;
}
}
return isThereNewIteration;
}
void TP11020::_checkInCodelets11021::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts11020[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iend1_darts11020[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->ist1_darts11020[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21i_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21im1_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31i_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31im1_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41i_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41im1_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51i_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51im1_darts11020[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts10796[this->getID()]);

/*printing node 11021: ForStmt*/
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
int** L2 = &(this->inputsTPParent->L2_darts11020[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts11020[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iend1 = &(this->inputsTPParent->iend1_darts11020[this->getLocalID()]);
(void)iend1/*OMP_SHARED_PRIVATE*/;
int** ist1 = &(this->inputsTPParent->ist1_darts11020[this->getLocalID()]);
(void)ist1/*OMP_SHARED_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts11020[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts11020[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts11020[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts11020[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21i = &(this->inputsTPParent->u21i_darts11020[this->getLocalID()]);
(void)u21i/*OMP_SHARED_PRIVATE*/;
double** u21im1 = &(this->inputsTPParent->u21im1_darts11020[this->getLocalID()]);
(void)u21im1/*OMP_SHARED_PRIVATE*/;
double** u31i = &(this->inputsTPParent->u31i_darts11020[this->getLocalID()]);
(void)u31i/*OMP_SHARED_PRIVATE*/;
double** u31im1 = &(this->inputsTPParent->u31im1_darts11020[this->getLocalID()]);
(void)u31im1/*OMP_SHARED_PRIVATE*/;
double** u41i = &(this->inputsTPParent->u41i_darts11020[this->getLocalID()]);
(void)u41i/*OMP_SHARED_PRIVATE*/;
double** u41im1 = &(this->inputsTPParent->u41im1_darts11020[this->getLocalID()]);
(void)u41im1/*OMP_SHARED_PRIVATE*/;
double** u51i = &(this->inputsTPParent->u51i_darts11020[this->getLocalID()]);
(void)u51i/*OMP_SHARED_PRIVATE*/;
double** u51im1 = &(this->inputsTPParent->u51im1_darts11020[this->getLocalID()]);
(void)u51im1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11020((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11020[0].decDep();
return;
}
for (int j_darts_counter_temp11020 = (*j);j_darts_counter_temp11020<=endRange && j_darts_counter_temp11020<=this->inputsTPParent->lastIteration11020;j_darts_counter_temp11020++)
{
{
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp11020 = (*k);
for(;k_darts_counter_temp11020 <= nz - 2;k_darts_counter_temp11020++){
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp11020 = (*i);
for(;i_darts_counter_temp11020 <= iend;i_darts_counter_temp11020++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11020 = (*m);
for(;m_darts_counter_temp11020 < 5;m_darts_counter_temp11020++){
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - tx2 * (flux[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - flux[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
}
(*m) = m_darts_counter_temp11020;
}
}
(*i) = i_darts_counter_temp11020;
}
(*(*L2)) = nx - 1;
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp11020 = (*i);
for(;i_darts_counter_temp11020 <= (*(*L2));i_darts_counter_temp11020++){
(*(*tmp)) = 1. / u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0];
(*(*u21i)) = (*(*tmp)) * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1];
(*(*u31i)) = (*(*tmp)) * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2];
(*(*u41i)) = (*(*tmp)) * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3];
(*(*u51i)) = (*(*tmp)) * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4];
(*(*tmp)) = 1. / u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0];
(*(*u21im1)) = (*(*tmp)) * u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1];
(*(*u31im1)) = (*(*tmp)) * u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2];
(*(*u41im1)) = (*(*tmp)) * u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3];
(*(*u51im1)) = (*(*tmp)) * u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4];
flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] = tx3 * ((*(*u31i)) - (*(*u31im1)));
flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] = tx3 * ((*(*u41i)) - (*(*u41im1)));
flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tx3 * (((((*(*u21i))) * ((*(*u21i)))) + (((*(*u31i))) * ((*(*u31i)))) + (((*(*u41i))) * ((*(*u41i))))) - ((((*(*u21im1))) * ((*(*u21im1)))) + (((*(*u31im1))) * ((*(*u31im1)))) + (((*(*u41im1))) * ((*(*u41im1)))))) + (1. / 6.) * tx3 * ((((*(*u21i))) * ((*(*u21i)))) - (((*(*u21im1))) * ((*(*u21im1))))) + 1.3999999999999999 * 1.3999999999999999 * tx3 * ((*(*u51i)) - (*(*u51im1)));
}
(*i) = i_darts_counter_temp11020;
}
{
/*Loop's init*/
(*i) = ist;
int i_darts_counter_temp11020 = (*i);
for(;i_darts_counter_temp11020 <= iend;i_darts_counter_temp11020++){
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0] + dx1 * tx1 * (u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0] - 2. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0] + u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][0]);
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] - flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1]) + dx2 * tx1 * (u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] - 2. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1] + u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][1]);
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] - flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2]) + dx3 * tx1 * (u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] - 2. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2] + u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][2]);
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] - flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3]) + dx4 * tx1 * (u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] - 2. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3] + u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][3]);
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] + tx3 * 0.10000000000000001 * 1. * (flux[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] - flux[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4]) + dx5 * tx1 * (u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] - 2. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4] + u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][4]);
}
(*i) = i_darts_counter_temp11020;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11020 = (*m);
for(;m_darts_counter_temp11020 < 5;m_darts_counter_temp11020++){
rsd[1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - dssp * (+5. * u[1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + u[3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
rsd[2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - dssp * (-4. * u[1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + 6. * u[2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + u[4][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
}
(*m) = m_darts_counter_temp11020;
}
(*(*ist1)) = 3;
(*(*iend1)) = nx - 4;
{
/*Loop's init*/
(*i) = (*(*ist1));
int i_darts_counter_temp11020 = (*i);
for(;i_darts_counter_temp11020 <= (*(*iend1));i_darts_counter_temp11020++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11020 = (*m);
for(;m_darts_counter_temp11020 < 5;m_darts_counter_temp11020++){
rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - dssp * (u[i_darts_counter_temp11020 - 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[i_darts_counter_temp11020 - 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + 6. * u[i_darts_counter_temp11020][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[i_darts_counter_temp11020 + 1][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + u[i_darts_counter_temp11020 + 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
}
(*m) = m_darts_counter_temp11020;
}
}
(*i) = i_darts_counter_temp11020;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11020 = (*m);
for(;m_darts_counter_temp11020 < 5;m_darts_counter_temp11020++){
rsd[nx - 3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[nx - 3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - dssp * (u[nx - 5][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[nx - 4][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + 6. * u[nx - 3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[nx - 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
rsd[nx - 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] = rsd[nx - 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - dssp * (u[nx - 4][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] - 4. * u[nx - 3][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020] + 5. * u[nx - 2][(j_darts_counter_temp11020)][k_darts_counter_temp11020][m_darts_counter_temp11020]);
}
(*m) = m_darts_counter_temp11020;
}
}
(*k) = k_darts_counter_temp11020;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11020[0].decDep();
}
TP11020::TP11020(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11020** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L2_darts11020(new int*[this->numThreads]),i_darts11020(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iend1_darts11020(new int*[this->numThreads]),ist1_darts11020(new int*[this->numThreads]),j_darts11020(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts11020(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts11020(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts11020(new double*[this->numThreads]),u21i_darts11020(new double*[this->numThreads]),u21im1_darts11020(new double*[this->numThreads]),u31i_darts11020(new double*[this->numThreads]),u31im1_darts11020(new double*[this->numThreads]),u41i_darts11020(new double*[this->numThreads]),u41im1_darts11020(new double*[this->numThreads]),u51i_darts11020(new double*[this->numThreads]),u51im1_darts11020(new double*[this->numThreads]), initIteration11020(in_initIteration), lastIteration11020(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets11020(new _barrierCodelets11020[1]) ,checkInCodelets11021(new _checkInCodelets11021[this->numThreads]){
/*Initialize the loop parameters*/
range11020 = abs (lastIteration11020 - initIteration11020) / 1;
rangePerCodelet11020 = range11020 / numThreads;
minIteration11020 = min<int>(lastIteration11020, initIteration11020);
remainderRange11020 = range11020 % numThreads;
/*Initialize inputs and vars.*/
this->L2_darts11020 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->iend1_darts11020 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->ist1_darts11020 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21i_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21im1_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31i_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31im1_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41i_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41im1_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51i_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51im1_darts11020 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets11020[0] = _barrierCodelets11020(this->numThreads,this->numThreads,this, 0);
_checkInCodelets11021 * checkInCodelets11021Ptr = (this->checkInCodelets11021);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets11021);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets11021Ptr) = _checkInCodelets11021(2,1,this,codeletCounter);
#else
(*checkInCodelets11021Ptr) = _checkInCodelets11021(1,1,this,codeletCounter);
#endif
(*checkInCodelets11021Ptr).decDep();
checkInCodelets11021Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP11020::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets11021[localID].setID (codeletID);
this->checkInCodelets11021[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets11021[localID + this->baseNumThreads * i] = _checkInCodelets11021(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets11021[localID + this->baseNumThreads * i] = _checkInCodelets11021(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets11021[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets11021[localID + this->baseNumThreads * i].decDep();
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
TP11020::~TP11020(){
delete [] L2_darts11020;
delete [] iend1_darts11020;
delete [] ist1_darts11020;
delete [] tmp_darts11020;
delete [] u21i_darts11020;
delete [] u21im1_darts11020;
delete [] u31i_darts11020;
delete [] u31im1_darts11020;
delete [] u41i_darts11020;
delete [] u41im1_darts11020;
delete [] u51i_darts11020;
delete [] u51im1_darts11020;
delete [] barrierCodelets11020;
delete [] checkInCodelets11021;
}
/*TP11660: OMPForDirective*/
void TP11660::_barrierCodelets11660::fire(void)
{
TP11660* myTP = static_cast<TP11660*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets11660[0].decDep ();
}
bool TP11660::requestNewRangeIterations11660(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet11660 * codeletID;
int tempEndRange   = rangePerCodelet11660 * (codeletID + 1);
if (remainderRange11660 != 0)
{
if (codeletID < (uint32_t)remainderRange11660)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange11660;
tempEndRange += remainderRange11660;
}
}
tempStartRange = tempStartRange*1 + minIteration11660;
tempEndRange = tempEndRange*1 + minIteration11660;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration11660 < lastIteration11660)
{
(this->inputsTPParent->i_darts11660[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts11660[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration11660;
}
}
return isThereNewIteration;
}
void TP11660::_checkInCodelets11661::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L1_darts11660[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts11660[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts11660[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31_darts11660[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31_darts10796[this->getID()]);

/*printing node 11661: ForStmt*/
/*var: L1*/
/*var: L2*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: q*/
/*var: u31*/
int** L1 = &(this->inputsTPParent->L1_darts11660[this->getLocalID()]);
(void)L1/*OMP_SHARED_PRIVATE*/;
int** L2 = &(this->inputsTPParent->L2_darts11660[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts11660[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts11660[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts11660[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts11660[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** u31 = &(this->inputsTPParent->u31_darts11660[this->getLocalID()]);
(void)u31/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11660((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11660[0].decDep();
return;
}
for (int i_darts_counter_temp11660 = (*i);i_darts_counter_temp11660<=endRange && i_darts_counter_temp11660<=this->inputsTPParent->lastIteration11660;i_darts_counter_temp11660++)
{
{
{
/*Loop's init*/
(*j) = (*(*L1));
int j_darts_counter_temp11660 = (*j);
for(;j_darts_counter_temp11660 <= (*(*L2));j_darts_counter_temp11660++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp11660 = (*k);
for(;k_darts_counter_temp11660 <= nz - 2;k_darts_counter_temp11660++){
flux[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][0] = u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2];
(*(*u31)) = u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2] / u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][0];
(*(*q)) = 0.5 * (u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][1] * u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][1] + u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2] * u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2] + u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][3] * u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][3]) / u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][0];
flux[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][1] = u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][1] * (*(*u31));
flux[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2] = u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][2] * (*(*u31)) + 0.40000000000000002 * (u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][4] - (*(*q)));
flux[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][3] = u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][3] * (*(*u31));
flux[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][4] = (1.3999999999999999 * u[(i_darts_counter_temp11660)][j_darts_counter_temp11660][k_darts_counter_temp11660][4] - 0.40000000000000002 * (*(*q))) * (*(*u31));
}
(*k) = k_darts_counter_temp11660;
}
}
(*j) = j_darts_counter_temp11660;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11660[0].decDep();
}
TP11660::TP11660(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11660** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L1_darts11660(new int*[this->numThreads]),L2_darts11660(new int*[this->numThreads]),i_darts11660(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts11660(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts11660(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts11660(new double*[this->numThreads]),u31_darts11660(new double*[this->numThreads]), initIteration11660(in_initIteration), lastIteration11660(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets11660(new _barrierCodelets11660[1]) ,checkInCodelets11661(new _checkInCodelets11661[this->numThreads]){
/*Initialize the loop parameters*/
range11660 = abs (lastIteration11660 - initIteration11660) / 1;
rangePerCodelet11660 = range11660 / numThreads;
minIteration11660 = min<int>(lastIteration11660, initIteration11660);
remainderRange11660 = range11660 % numThreads;
/*Initialize inputs and vars.*/
this->L1_darts11660 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->L2_darts11660 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->q_darts11660 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31_darts11660 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets11660[0] = _barrierCodelets11660(this->numThreads,this->numThreads,this, 0);
_checkInCodelets11661 * checkInCodelets11661Ptr = (this->checkInCodelets11661);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets11661);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets11661Ptr) = _checkInCodelets11661(2,1,this,codeletCounter);
#else
(*checkInCodelets11661Ptr) = _checkInCodelets11661(1,1,this,codeletCounter);
#endif
(*checkInCodelets11661Ptr).decDep();
checkInCodelets11661Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP11660::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets11661[localID].setID (codeletID);
this->checkInCodelets11661[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets11661[localID + this->baseNumThreads * i] = _checkInCodelets11661(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets11661[localID + this->baseNumThreads * i] = _checkInCodelets11661(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets11661[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets11661[localID + this->baseNumThreads * i].decDep();
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
TP11660::~TP11660(){
delete [] L1_darts11660;
delete [] L2_darts11660;
delete [] q_darts11660;
delete [] u31_darts11660;
delete [] barrierCodelets11660;
delete [] checkInCodelets11661;
}
/*TP11809: OMPForDirective*/
void TP11809::_barrierCodelets11809::fire(void)
{
TP11809* myTP = static_cast<TP11809*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets11809[0].decDep ();
}
bool TP11809::requestNewRangeIterations11809(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet11809 * codeletID;
int tempEndRange   = rangePerCodelet11809 * (codeletID + 1);
if (remainderRange11809 != 0)
{
if (codeletID < (uint32_t)remainderRange11809)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange11809;
tempEndRange += remainderRange11809;
}
}
tempStartRange = tempStartRange*1 + minIteration11809;
tempEndRange = tempEndRange*1 + minIteration11809;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration11809 < lastIteration11809)
{
(this->inputsTPParent->i_darts11809[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts11809[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration11809;
}
}
return isThereNewIteration;
}
void TP11809::_checkInCodelets11810::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->L2_darts11809[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jend1_darts11809[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jst1_darts11809[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21j_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21jm1_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31j_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31jm1_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41j_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41jm1_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51j_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51jm1_darts11809[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts10796[this->getID()]);

/*printing node 11810: ForStmt*/
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
int** L2 = &(this->inputsTPParent->L2_darts11809[this->getLocalID()]);
(void)L2/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts11809[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts11809[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jend1 = &(this->inputsTPParent->jend1_darts11809[this->getLocalID()]);
(void)jend1/*OMP_SHARED_PRIVATE*/;
int** jst1 = &(this->inputsTPParent->jst1_darts11809[this->getLocalID()]);
(void)jst1/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts11809[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts11809[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts11809[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21j = &(this->inputsTPParent->u21j_darts11809[this->getLocalID()]);
(void)u21j/*OMP_SHARED_PRIVATE*/;
double** u21jm1 = &(this->inputsTPParent->u21jm1_darts11809[this->getLocalID()]);
(void)u21jm1/*OMP_SHARED_PRIVATE*/;
double** u31j = &(this->inputsTPParent->u31j_darts11809[this->getLocalID()]);
(void)u31j/*OMP_SHARED_PRIVATE*/;
double** u31jm1 = &(this->inputsTPParent->u31jm1_darts11809[this->getLocalID()]);
(void)u31jm1/*OMP_SHARED_PRIVATE*/;
double** u41j = &(this->inputsTPParent->u41j_darts11809[this->getLocalID()]);
(void)u41j/*OMP_SHARED_PRIVATE*/;
double** u41jm1 = &(this->inputsTPParent->u41jm1_darts11809[this->getLocalID()]);
(void)u41jm1/*OMP_SHARED_PRIVATE*/;
double** u51j = &(this->inputsTPParent->u51j_darts11809[this->getLocalID()]);
(void)u51j/*OMP_SHARED_PRIVATE*/;
double** u51jm1 = &(this->inputsTPParent->u51jm1_darts11809[this->getLocalID()]);
(void)u51jm1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11809((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11809[0].decDep();
return;
}
for (int i_darts_counter_temp11809 = (*i);i_darts_counter_temp11809<=endRange && i_darts_counter_temp11809<=this->inputsTPParent->lastIteration11809;i_darts_counter_temp11809++)
{
{
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp11809 = (*k);
for(;k_darts_counter_temp11809 <= nz - 2;k_darts_counter_temp11809++){
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp11809 = (*j);
for(;j_darts_counter_temp11809 <= jend;j_darts_counter_temp11809++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11809 = (*m);
for(;m_darts_counter_temp11809 < 5;m_darts_counter_temp11809++){
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][m_darts_counter_temp11809] - ty2 * (flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][m_darts_counter_temp11809] - flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][m_darts_counter_temp11809]);
}
(*m) = m_darts_counter_temp11809;
}
}
(*j) = j_darts_counter_temp11809;
}
(*(*L2)) = ny - 1;
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp11809 = (*j);
for(;j_darts_counter_temp11809 <= (*(*L2));j_darts_counter_temp11809++){
(*(*tmp)) = 1. / u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][0];
(*(*u21j)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1];
(*(*u31j)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2];
(*(*u41j)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3];
(*(*u51j)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4];
(*(*tmp)) = 1. / u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][0];
(*(*u21jm1)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][1];
(*(*u31jm1)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][2];
(*(*u41jm1)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][3];
(*(*u51jm1)) = (*(*tmp)) * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][4];
flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1] = ty3 * ((*(*u21j)) - (*(*u21jm1)));
flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2] = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3] = ty3 * ((*(*u41j)) - (*(*u41jm1)));
flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * ty3 * (((((*(*u21j))) * ((*(*u21j)))) + (((*(*u31j))) * ((*(*u31j)))) + (((*(*u41j))) * ((*(*u41j))))) - ((((*(*u21jm1))) * ((*(*u21jm1)))) + (((*(*u31jm1))) * ((*(*u31jm1)))) + (((*(*u41jm1))) * ((*(*u41jm1)))))) + (1. / 6.) * ty3 * ((((*(*u31j))) * ((*(*u31j)))) - (((*(*u31jm1))) * ((*(*u31jm1))))) + 1.3999999999999999 * 1.3999999999999999 * ty3 * ((*(*u51j)) - (*(*u51jm1)));
}
(*j) = j_darts_counter_temp11809;
}
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp11809 = (*j);
for(;j_darts_counter_temp11809 <= jend;j_darts_counter_temp11809++){
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][0] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][0] + dy1 * ty1 * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][0] - 2. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][0] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][0]);
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][1] - flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1]) + dy2 * ty1 * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][1] - 2. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][1] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][1]);
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][2] - flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2]) + dy3 * ty1 * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][2] - 2. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][2] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][2]);
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][3] - flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3]) + dy4 * ty1 * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][3] - 2. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][3] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][3]);
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4] + ty3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][4] - flux[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4]) + dy5 * ty1 * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][4] - 2. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][4] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][4]);
}
(*j) = j_darts_counter_temp11809;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11809 = (*m);
for(;m_darts_counter_temp11809 < 5;m_darts_counter_temp11809++){
rsd[(i_darts_counter_temp11809)][1][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][1][k_darts_counter_temp11809][m_darts_counter_temp11809] - dssp * (+5. * u[(i_darts_counter_temp11809)][1][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][2][k_darts_counter_temp11809][m_darts_counter_temp11809] + u[(i_darts_counter_temp11809)][3][k_darts_counter_temp11809][m_darts_counter_temp11809]);
rsd[(i_darts_counter_temp11809)][2][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][2][k_darts_counter_temp11809][m_darts_counter_temp11809] - dssp * (-4. * u[(i_darts_counter_temp11809)][1][k_darts_counter_temp11809][m_darts_counter_temp11809] + 6. * u[(i_darts_counter_temp11809)][2][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][3][k_darts_counter_temp11809][m_darts_counter_temp11809] + u[(i_darts_counter_temp11809)][4][k_darts_counter_temp11809][m_darts_counter_temp11809]);
}
(*m) = m_darts_counter_temp11809;
}
(*(*jst1)) = 3;
(*(*jend1)) = ny - 4;
{
/*Loop's init*/
(*j) = (*(*jst1));
int j_darts_counter_temp11809 = (*j);
for(;j_darts_counter_temp11809 <= (*(*jend1));j_darts_counter_temp11809++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11809 = (*m);
for(;m_darts_counter_temp11809 < 5;m_darts_counter_temp11809++){
rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][m_darts_counter_temp11809] - dssp * (u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 2][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 - 1][k_darts_counter_temp11809][m_darts_counter_temp11809] + 6. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 1][k_darts_counter_temp11809][m_darts_counter_temp11809] + u[(i_darts_counter_temp11809)][j_darts_counter_temp11809 + 2][k_darts_counter_temp11809][m_darts_counter_temp11809]);
}
(*m) = m_darts_counter_temp11809;
}
}
(*j) = j_darts_counter_temp11809;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp11809 = (*m);
for(;m_darts_counter_temp11809 < 5;m_darts_counter_temp11809++){
rsd[(i_darts_counter_temp11809)][ny - 3][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][ny - 3][k_darts_counter_temp11809][m_darts_counter_temp11809] - dssp * (u[(i_darts_counter_temp11809)][ny - 5][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][ny - 4][k_darts_counter_temp11809][m_darts_counter_temp11809] + 6. * u[(i_darts_counter_temp11809)][ny - 3][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][ny - 2][k_darts_counter_temp11809][m_darts_counter_temp11809]);
rsd[(i_darts_counter_temp11809)][ny - 2][k_darts_counter_temp11809][m_darts_counter_temp11809] = rsd[(i_darts_counter_temp11809)][ny - 2][k_darts_counter_temp11809][m_darts_counter_temp11809] - dssp * (u[(i_darts_counter_temp11809)][ny - 4][k_darts_counter_temp11809][m_darts_counter_temp11809] - 4. * u[(i_darts_counter_temp11809)][ny - 3][k_darts_counter_temp11809][m_darts_counter_temp11809] + 5. * u[(i_darts_counter_temp11809)][ny - 2][k_darts_counter_temp11809][m_darts_counter_temp11809]);
}
(*m) = m_darts_counter_temp11809;
}
}
(*k) = k_darts_counter_temp11809;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets11809[0].decDep();
}
TP11809::TP11809(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11809** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),L2_darts11809(new int*[this->numThreads]),i_darts11809(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts11809(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jend1_darts11809(new int*[this->numThreads]),jst1_darts11809(new int*[this->numThreads]),k_darts11809(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts11809(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts11809(new double*[this->numThreads]),u21j_darts11809(new double*[this->numThreads]),u21jm1_darts11809(new double*[this->numThreads]),u31j_darts11809(new double*[this->numThreads]),u31jm1_darts11809(new double*[this->numThreads]),u41j_darts11809(new double*[this->numThreads]),u41jm1_darts11809(new double*[this->numThreads]),u51j_darts11809(new double*[this->numThreads]),u51jm1_darts11809(new double*[this->numThreads]), initIteration11809(in_initIteration), lastIteration11809(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets11809(new _barrierCodelets11809[1]) ,checkInCodelets11810(new _checkInCodelets11810[this->numThreads]){
/*Initialize the loop parameters*/
range11809 = abs (lastIteration11809 - initIteration11809) / 1;
rangePerCodelet11809 = range11809 / numThreads;
minIteration11809 = min<int>(lastIteration11809, initIteration11809);
remainderRange11809 = range11809 % numThreads;
/*Initialize inputs and vars.*/
this->L2_darts11809 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jend1_darts11809 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jst1_darts11809 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21j_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21jm1_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31j_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31jm1_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41j_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41jm1_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51j_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51jm1_darts11809 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets11809[0] = _barrierCodelets11809(this->numThreads,this->numThreads,this, 0);
_checkInCodelets11810 * checkInCodelets11810Ptr = (this->checkInCodelets11810);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets11810);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets11810Ptr) = _checkInCodelets11810(2,1,this,codeletCounter);
#else
(*checkInCodelets11810Ptr) = _checkInCodelets11810(1,1,this,codeletCounter);
#endif
(*checkInCodelets11810Ptr).decDep();
checkInCodelets11810Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP11809::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets11810[localID].setID (codeletID);
this->checkInCodelets11810[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets11810[localID + this->baseNumThreads * i] = _checkInCodelets11810(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets11810[localID + this->baseNumThreads * i] = _checkInCodelets11810(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets11810[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets11810[localID + this->baseNumThreads * i].decDep();
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
TP11809::~TP11809(){
delete [] L2_darts11809;
delete [] jend1_darts11809;
delete [] jst1_darts11809;
delete [] tmp_darts11809;
delete [] u21j_darts11809;
delete [] u21jm1_darts11809;
delete [] u31j_darts11809;
delete [] u31jm1_darts11809;
delete [] u41j_darts11809;
delete [] u41jm1_darts11809;
delete [] u51j_darts11809;
delete [] u51jm1_darts11809;
delete [] barrierCodelets11809;
delete [] checkInCodelets11810;
}
/*TP12446: OMPForDirective*/
void TP12446::_barrierCodelets12446::fire(void)
{
TP12446* myTP = static_cast<TP12446*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets12446[0].decDep ();
}
bool TP12446::requestNewRangeIterations12446(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet12446 * codeletID;
int tempEndRange   = rangePerCodelet12446 * (codeletID + 1);
if (remainderRange12446 != 0)
{
if (codeletID < (uint32_t)remainderRange12446)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange12446;
tempEndRange += remainderRange12446;
}
}
tempStartRange = tempStartRange*1 + minIteration12446;
tempEndRange = tempEndRange*1 + minIteration12446;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration12446 < lastIteration12446)
{
(this->inputsTPParent->i_darts12446[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts12446[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration12446;
}
}
return isThereNewIteration;
}
void TP12446::_checkInCodelets12447::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->q_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->q_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->tmp_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21k_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u21km1_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31k_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u31km1_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41k_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u41km1_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51k_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts10796[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->u51km1_darts12446[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts10796[this->getID()]);

/*printing node 12447: ForStmt*/
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
int* i = &(this->inputsTPParent->i_darts12446[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts12446[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts12446[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts12446[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** q = &(this->inputsTPParent->q_darts12446[this->getLocalID()]);
(void)q/*OMP_SHARED_PRIVATE*/;
double** tmp = &(this->inputsTPParent->tmp_darts12446[this->getLocalID()]);
(void)tmp/*OMP_SHARED_PRIVATE*/;
double** u21k = &(this->inputsTPParent->u21k_darts12446[this->getLocalID()]);
(void)u21k/*OMP_SHARED_PRIVATE*/;
double** u21km1 = &(this->inputsTPParent->u21km1_darts12446[this->getLocalID()]);
(void)u21km1/*OMP_SHARED_PRIVATE*/;
double** u31k = &(this->inputsTPParent->u31k_darts12446[this->getLocalID()]);
(void)u31k/*OMP_SHARED_PRIVATE*/;
double** u31km1 = &(this->inputsTPParent->u31km1_darts12446[this->getLocalID()]);
(void)u31km1/*OMP_SHARED_PRIVATE*/;
double** u41 = &(this->inputsTPParent->u41_darts12446[this->getLocalID()]);
(void)u41/*OMP_SHARED_PRIVATE*/;
double** u41k = &(this->inputsTPParent->u41k_darts12446[this->getLocalID()]);
(void)u41k/*OMP_SHARED_PRIVATE*/;
double** u41km1 = &(this->inputsTPParent->u41km1_darts12446[this->getLocalID()]);
(void)u41km1/*OMP_SHARED_PRIVATE*/;
double** u51k = &(this->inputsTPParent->u51k_darts12446[this->getLocalID()]);
(void)u51k/*OMP_SHARED_PRIVATE*/;
double** u51km1 = &(this->inputsTPParent->u51km1_darts12446[this->getLocalID()]);
(void)u51km1/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations12446((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets12446[0].decDep();
return;
}
for (int i_darts_counter_temp12446 = (*i);i_darts_counter_temp12446<=endRange && i_darts_counter_temp12446<=this->inputsTPParent->lastIteration12446;i_darts_counter_temp12446++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp12446 = (*j);
for(;j_darts_counter_temp12446 <= jend;j_darts_counter_temp12446++){
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp12446 = (*k);
for(;k_darts_counter_temp12446 <= nz - 1;k_darts_counter_temp12446++){
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0] = u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3];
(*(*u41)) = u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] / u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0];
(*(*q)) = 0.5 * (u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3]) / u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0];
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] = u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] * (*(*u41));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] = u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] * (*(*u41));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] = u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] * (*(*u41)) + 0.40000000000000002 * (u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] - (*(*q)));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] = (1.3999999999999999 * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] - 0.40000000000000002 * (*(*q))) * (*(*u41));
}
(*k) = k_darts_counter_temp12446;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp12446 = (*k);
for(;k_darts_counter_temp12446 <= nz - 2;k_darts_counter_temp12446++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp12446 = (*m);
for(;m_darts_counter_temp12446 < 5;m_darts_counter_temp12446++){
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][m_darts_counter_temp12446] - tz2 * (flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][m_darts_counter_temp12446] - flux[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][m_darts_counter_temp12446]);
}
(*m) = m_darts_counter_temp12446;
}
}
(*k) = k_darts_counter_temp12446;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp12446 = (*k);
for(;k_darts_counter_temp12446 <= nz - 1;k_darts_counter_temp12446++){
(*(*tmp)) = 1. / u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0];
(*(*u21k)) = (*(*tmp)) * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1];
(*(*u31k)) = (*(*tmp)) * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2];
(*(*u41k)) = (*(*tmp)) * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3];
(*(*u51k)) = (*(*tmp)) * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4];
(*(*tmp)) = 1. / u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][0];
(*(*u21km1)) = (*(*tmp)) * u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][1];
(*(*u31km1)) = (*(*tmp)) * u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][2];
(*(*u41km1)) = (*(*tmp)) * u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][3];
(*(*u51km1)) = (*(*tmp)) * u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][4];
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] = tz3 * ((*(*u21k)) - (*(*u21km1)));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] = tz3 * ((*(*u31k)) - (*(*u31km1)));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tz3 * (((((*(*u21k))) * ((*(*u21k)))) + (((*(*u31k))) * ((*(*u31k)))) + (((*(*u41k))) * ((*(*u41k))))) - ((((*(*u21km1))) * ((*(*u21km1)))) + (((*(*u31km1))) * ((*(*u31km1)))) + (((*(*u41km1))) * ((*(*u41km1)))))) + (1. / 6.) * tz3 * ((((*(*u41k))) * ((*(*u41k)))) - (((*(*u41km1))) * ((*(*u41km1))))) + 1.3999999999999999 * 1.3999999999999999 * tz3 * ((*(*u51k)) - (*(*u51km1)));
}
(*k) = k_darts_counter_temp12446;
}
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp12446 = (*k);
for(;k_darts_counter_temp12446 <= nz - 2;k_darts_counter_temp12446++){
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0] + dz1 * tz1 * (u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][0] - 2. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][0] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][0]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][1] - flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1]) + dz2 * tz1 * (u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][1] - 2. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][1] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][1]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][2] - flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2]) + dz3 * tz1 * (u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][2] - 2. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][2] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][2]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][3] - flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3]) + dz4 * tz1 * (u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][3] - 2. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][3] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][3]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] + tz3 * 0.10000000000000001 * 1. * (flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][4] - flux[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4]) + dz5 * tz1 * (u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][4] - 2. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][4] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][4]);
}
(*k) = k_darts_counter_temp12446;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp12446 = (*m);
for(;m_darts_counter_temp12446 < 5;m_darts_counter_temp12446++){
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][1][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][1][m_darts_counter_temp12446] - dssp * (+5. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][1][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][2][m_darts_counter_temp12446] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][3][m_darts_counter_temp12446]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][2][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][2][m_darts_counter_temp12446] - dssp * (-4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][1][m_darts_counter_temp12446] + 6. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][2][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][3][m_darts_counter_temp12446] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][4][m_darts_counter_temp12446]);
}
(*m) = m_darts_counter_temp12446;
}
{
/*Loop's init*/
(*k) = 3;
int k_darts_counter_temp12446 = (*k);
for(;k_darts_counter_temp12446 <= nz - 4;k_darts_counter_temp12446++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp12446 = (*m);
for(;m_darts_counter_temp12446 < 5;m_darts_counter_temp12446++){
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][m_darts_counter_temp12446] - dssp * (u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 - 2][m_darts_counter_temp12446] - 4. * u[k_darts_counter_temp12446 - 1][(i_darts_counter_temp12446)][j_darts_counter_temp12446][m_darts_counter_temp12446] + 6. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 1][m_darts_counter_temp12446] + u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][k_darts_counter_temp12446 + 2][m_darts_counter_temp12446]);
}
(*m) = m_darts_counter_temp12446;
}
}
(*k) = k_darts_counter_temp12446;
}
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp12446 = (*m);
for(;m_darts_counter_temp12446 < 5;m_darts_counter_temp12446++){
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 3][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 3][m_darts_counter_temp12446] - dssp * (u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 5][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 4][m_darts_counter_temp12446] + 6. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 3][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 2][m_darts_counter_temp12446]);
rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 2][m_darts_counter_temp12446] = rsd[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 2][m_darts_counter_temp12446] - dssp * (u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 4][m_darts_counter_temp12446] - 4. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 3][m_darts_counter_temp12446] + 5. * u[(i_darts_counter_temp12446)][j_darts_counter_temp12446][nz - 2][m_darts_counter_temp12446]);
}
(*m) = m_darts_counter_temp12446;
}
}
(*j) = j_darts_counter_temp12446;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets12446[0].decDep();
}
TP12446::TP12446(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP12446** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts12446(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts12446(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts12446(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts12446(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,q_darts12446(new double*[this->numThreads]),tmp_darts12446(new double*[this->numThreads]),u21k_darts12446(new double*[this->numThreads]),u21km1_darts12446(new double*[this->numThreads]),u31k_darts12446(new double*[this->numThreads]),u31km1_darts12446(new double*[this->numThreads]),u41_darts12446(new double*[this->numThreads]),u41k_darts12446(new double*[this->numThreads]),u41km1_darts12446(new double*[this->numThreads]),u51k_darts12446(new double*[this->numThreads]),u51km1_darts12446(new double*[this->numThreads]), initIteration12446(in_initIteration), lastIteration12446(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets12446(new _barrierCodelets12446[1]) ,checkInCodelets12447(new _checkInCodelets12447[this->numThreads]){
/*Initialize the loop parameters*/
range12446 = abs (lastIteration12446 - initIteration12446) / 1;
rangePerCodelet12446 = range12446 / numThreads;
minIteration12446 = min<int>(lastIteration12446, initIteration12446);
remainderRange12446 = range12446 % numThreads;
/*Initialize inputs and vars.*/
this->q_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->tmp_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21k_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u21km1_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31k_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u31km1_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41k_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u41km1_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51k_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->u51km1_darts12446 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets12446[0] = _barrierCodelets12446(this->numThreads,this->numThreads,this, 0);
_checkInCodelets12447 * checkInCodelets12447Ptr = (this->checkInCodelets12447);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets12447);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets12447Ptr) = _checkInCodelets12447(2,1,this,codeletCounter);
#else
(*checkInCodelets12447Ptr) = _checkInCodelets12447(1,1,this,codeletCounter);
#endif
(*checkInCodelets12447Ptr).decDep();
checkInCodelets12447Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP12446::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets12447[localID].setID (codeletID);
this->checkInCodelets12447[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets12447[localID + this->baseNumThreads * i] = _checkInCodelets12447(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets12447[localID + this->baseNumThreads * i] = _checkInCodelets12447(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets12447[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets12447[localID + this->baseNumThreads * i].decDep();
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
TP12446::~TP12446(){
delete [] q_darts12446;
delete [] tmp_darts12446;
delete [] u21k_darts12446;
delete [] u21km1_darts12446;
delete [] u31k_darts12446;
delete [] u31km1_darts12446;
delete [] u41_darts12446;
delete [] u41k_darts12446;
delete [] u41km1_darts12446;
delete [] u51k_darts12446;
delete [] u51km1_darts12446;
delete [] barrierCodelets12446;
delete [] checkInCodelets12447;
}
/*TP13197: OMPParallelDirective*/
void TP13197::_barrierCodelets13197::fire(void)
{
TP13197* myTP = static_cast<TP13197*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP13197::_checkInCodelets13201::fire(void)
{
/*region 13201 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13201;
if(idx < myTP->TPsToUse13201){
if (!__sync_val_compare_and_swap (&(myTP->TP13201_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13201;
int minIteration = min<int >(nx, 0);
int remainderRange = range % myTP->TPsToUse13201;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13201 - 1)
{
lastIteration = nx;
}
#if USEINVOKE == 1
invoke < TP13201 > (myTP, myTP->codeletsPerTP13201 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13201Ptr[idx]));
#else
place < TP13201 > (idx, myTP, myTP->codeletsPerTP13201 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13201Ptr[idx]));
#endif
}else{
if (myTP->TP13201Ptr[idx] != nullptr){
myTP->TP13201Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13197::_barrierCodelets13201::fire(void)
{
TP13197* myTP =  static_cast<TP13197*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets13252[codeletsCounter].decDep();
}
}
}
void TP13197::_checkInCodelets13252::fire(void)
{
/*region 13252 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13252;
if(idx < myTP->TPsToUse13252){
if (!__sync_val_compare_and_swap (&(myTP->TP13252_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13252;
int minIteration = min<int >(nx, 0);
int remainderRange = range % myTP->TPsToUse13252;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13252 - 1)
{
lastIteration = nx;
}
#if USEINVOKE == 1
invoke < TP13252 > (myTP, myTP->codeletsPerTP13252 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13252Ptr[idx]));
#else
place < TP13252 > (idx, myTP, myTP->codeletsPerTP13252 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13252Ptr[idx]));
#endif
}else{
if (myTP->TP13252Ptr[idx] != nullptr){
myTP->TP13252Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13197::_barrierCodelets13252::fire(void)
{
TP13197* myTP =  static_cast<TP13197*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets13294[codeletsCounter].decDep();
}
}
}
void TP13197::_checkInCodelets13294::fire(void)
{
/*region 13294 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13294;
if(idx < myTP->TPsToUse13294){
if (!__sync_val_compare_and_swap (&(myTP->TP13294_alreadyLaunched[idx]), 0, 1)){
int range = abs (nx - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13294;
int minIteration = min<int >(nx, 0);
int remainderRange = range % myTP->TPsToUse13294;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < nx)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13294 - 1)
{
lastIteration = nx;
}
#if USEINVOKE == 1
invoke < TP13294 > (myTP, myTP->codeletsPerTP13294 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13294Ptr[idx]));
#else
place < TP13294 > (idx, myTP, myTP->codeletsPerTP13294 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13294Ptr[idx]));
#endif
}else{
if (myTP->TP13294Ptr[idx] != nullptr){
myTP->TP13294Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13197::_barrierCodelets13294::fire(void)
{
TP13197* myTP =  static_cast<TP13197*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets13338[codeletsCounter].decDep();
}
}
}
void TP13197::_checkInCodelets13338::fire(void)
{
/*region 13338 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13338;
if(idx < myTP->TPsToUse13338){
if (!__sync_val_compare_and_swap (&(myTP->TP13338_alreadyLaunched[idx]), 0, 1)){
int range = abs (ny - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13338;
int minIteration = min<int >(ny, 0);
int remainderRange = range % myTP->TPsToUse13338;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < ny)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13338 - 1)
{
lastIteration = ny;
}
#if USEINVOKE == 1
invoke < TP13338 > (myTP, myTP->codeletsPerTP13338 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13338Ptr[idx]));
#else
place < TP13338 > (idx, myTP, myTP->codeletsPerTP13338 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13338Ptr[idx]));
#endif
}else{
if (myTP->TP13338Ptr[idx] != nullptr){
myTP->TP13338Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13197::_barrierCodelets13338::fire(void)
{
TP13197* myTP =  static_cast<TP13197*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets13380[codeletsCounter].decDep();
}
}
}
void TP13197::_checkInCodelets13380::fire(void)
{
/*region 13380 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13380;
if(idx < myTP->TPsToUse13380){
if (!__sync_val_compare_and_swap (&(myTP->TP13380_alreadyLaunched[idx]), 0, 1)){
int range = abs (ny - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13380;
int minIteration = min<int >(ny, 0);
int remainderRange = range % myTP->TPsToUse13380;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < ny)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13380 - 1)
{
lastIteration = ny;
}
#if USEINVOKE == 1
invoke < TP13380 > (myTP, myTP->codeletsPerTP13380 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13380Ptr[idx]));
#else
place < TP13380 > (idx, myTP, myTP->codeletsPerTP13380 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13380Ptr[idx]));
#endif
}else{
if (myTP->TP13380Ptr[idx] != nullptr){
myTP->TP13380Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13197::_barrierCodelets13380::fire(void)
{
TP13197* myTP =  static_cast<TP13197*>(myTP_);
myTP->TPParent->barrierCodelets13197[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets13197[0]));
}
TP13197::TP13197(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),i_darts13197(new int[this->numThreads])/*VARIABLE*/,iglob_darts13197(new int[this->numThreads])/*VARIABLE*/,j_darts13197(new int[this->numThreads])/*VARIABLE*/,jglob_darts13197(new int[this->numThreads])/*VARIABLE*/,k_darts13197(new int[this->numThreads])/*VARIABLE*/, TP13201Ptr(new TP13201 *[NUMTPS13201]), TP13201_alreadyLaunched(new size_t [NUMTPS13201]), numTPsSet13201(0), numTPsReady13201(0), TPsToUse13201(NUMTPS13201), codeletsPerTP13201(this->numThreads/NUMTPS13201), totalCodelets13201(this->TPsToUse13201*this->codeletsPerTP13201), TP13252Ptr(new TP13252 *[NUMTPS13252]), TP13252_alreadyLaunched(new size_t [NUMTPS13252]), numTPsSet13252(0), numTPsReady13252(0), TPsToUse13252(NUMTPS13252), codeletsPerTP13252(this->numThreads/NUMTPS13252), totalCodelets13252(this->TPsToUse13252*this->codeletsPerTP13252), TP13294Ptr(new TP13294 *[NUMTPS13294]), TP13294_alreadyLaunched(new size_t [NUMTPS13294]), numTPsSet13294(0), numTPsReady13294(0), TPsToUse13294(NUMTPS13294), codeletsPerTP13294(this->numThreads/NUMTPS13294), totalCodelets13294(this->TPsToUse13294*this->codeletsPerTP13294), TP13338Ptr(new TP13338 *[NUMTPS13338]), TP13338_alreadyLaunched(new size_t [NUMTPS13338]), numTPsSet13338(0), numTPsReady13338(0), TPsToUse13338(NUMTPS13338), codeletsPerTP13338(this->numThreads/NUMTPS13338), totalCodelets13338(this->TPsToUse13338*this->codeletsPerTP13338), TP13380Ptr(new TP13380 *[NUMTPS13380]), TP13380_alreadyLaunched(new size_t [NUMTPS13380]), numTPsSet13380(0), numTPsReady13380(0), TPsToUse13380(NUMTPS13380), codeletsPerTP13380(this->numThreads/NUMTPS13380), totalCodelets13380(this->TPsToUse13380*this->codeletsPerTP13380) ,barrierCodelets13197(new _barrierCodelets13197[1]) ,checkInCodelets13201(new _checkInCodelets13201[this->numThreads]) ,barrierCodelets13201(new _barrierCodelets13201[1]) ,checkInCodelets13252(new _checkInCodelets13252[this->numThreads]) ,barrierCodelets13252(new _barrierCodelets13252[1]) ,checkInCodelets13294(new _checkInCodelets13294[this->numThreads]) ,barrierCodelets13294(new _barrierCodelets13294[1]) ,checkInCodelets13338(new _checkInCodelets13338[this->numThreads]) ,barrierCodelets13338(new _barrierCodelets13338[1]) ,checkInCodelets13380(new _checkInCodelets13380[this->numThreads]) ,barrierCodelets13380(new _barrierCodelets13380[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13197[0] = _barrierCodelets13197(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets13380[0] = _barrierCodelets13380(NUMTPS13380,NUMTPS13380,this, 0);
barrierCodelets13338[0] = _barrierCodelets13338(NUMTPS13338,NUMTPS13338,this, 0);
barrierCodelets13294[0] = _barrierCodelets13294(NUMTPS13294,NUMTPS13294,this, 0);
barrierCodelets13252[0] = _barrierCodelets13252(NUMTPS13252,NUMTPS13252,this, 0);
barrierCodelets13201[0] = _barrierCodelets13201(NUMTPS13201,NUMTPS13201,this, 0);
_checkInCodelets13380 * checkInCodelets13380Ptr = (this->checkInCodelets13380);
for(int i=0; i<NUMTPS13380; i++)
{
TP13380Ptr[i] = nullptr;
TP13380_alreadyLaunched[i] = 0;
}
_checkInCodelets13338 * checkInCodelets13338Ptr = (this->checkInCodelets13338);
for(int i=0; i<NUMTPS13338; i++)
{
TP13338Ptr[i] = nullptr;
TP13338_alreadyLaunched[i] = 0;
}
_checkInCodelets13294 * checkInCodelets13294Ptr = (this->checkInCodelets13294);
for(int i=0; i<NUMTPS13294; i++)
{
TP13294Ptr[i] = nullptr;
TP13294_alreadyLaunched[i] = 0;
}
_checkInCodelets13252 * checkInCodelets13252Ptr = (this->checkInCodelets13252);
for(int i=0; i<NUMTPS13252; i++)
{
TP13252Ptr[i] = nullptr;
TP13252_alreadyLaunched[i] = 0;
}
_checkInCodelets13201 * checkInCodelets13201Ptr = (this->checkInCodelets13201);
for(int i=0; i<NUMTPS13201; i++)
{
TP13201Ptr[i] = nullptr;
TP13201_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets13380Ptr) = _checkInCodelets13380(1,1,this,codeletCounter);
checkInCodelets13380Ptr++;
(*checkInCodelets13338Ptr) = _checkInCodelets13338(1,1,this,codeletCounter);
checkInCodelets13338Ptr++;
(*checkInCodelets13294Ptr) = _checkInCodelets13294(1,1,this,codeletCounter);
checkInCodelets13294Ptr++;
(*checkInCodelets13252Ptr) = _checkInCodelets13252(1,1,this,codeletCounter);
checkInCodelets13252Ptr++;
(*checkInCodelets13201Ptr) = _checkInCodelets13201(1,1,this,codeletCounter);
(*checkInCodelets13201Ptr).decDep();
checkInCodelets13201Ptr++;
}
}
TP13197::~TP13197(){
delete []i_darts13197;
delete []iglob_darts13197;
delete []j_darts13197;
delete []jglob_darts13197;
delete []k_darts13197;
delete [] barrierCodelets13197;
delete [] barrierCodelets13380;
delete [] checkInCodelets13380;
delete [] barrierCodelets13338;
delete [] checkInCodelets13338;
delete [] barrierCodelets13294;
delete [] checkInCodelets13294;
delete [] barrierCodelets13252;
delete [] checkInCodelets13252;
delete [] barrierCodelets13201;
delete [] checkInCodelets13201;
}
/*TP13201: OMPForDirective*/
void TP13201::_barrierCodelets13201::fire(void)
{
TP13201* myTP = static_cast<TP13201*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13201[0].decDep ();
}
bool TP13201::requestNewRangeIterations13201(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13201 * codeletID;
int tempEndRange   = rangePerCodelet13201 * (codeletID + 1);
if (remainderRange13201 != 0)
{
if (codeletID < (uint32_t)remainderRange13201)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13201;
tempEndRange += remainderRange13201;
}
}
tempStartRange = tempStartRange*1 + minIteration13201;
tempEndRange = tempEndRange*1 + minIteration13201;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13201 < lastIteration13201)
{
(this->inputsTPParent->i_darts13201[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts13201[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13201;
}
}
return isThereNewIteration;
}
void TP13201::_checkInCodelets13202::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iglob_darts13201[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13197[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jglob_darts13201[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13197[this->getID()]);

/*printing node 13202: ForStmt*/
/*var: i*/
/*var: iglob*/
/*var: j*/
/*var: jglob*/
int* i = &(this->inputsTPParent->i_darts13201[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iglob = &(this->inputsTPParent->iglob_darts13201[this->getLocalID()]);
(void)iglob/*OMP_SHARED_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts13201[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jglob = &(this->inputsTPParent->jglob_darts13201[this->getLocalID()]);
(void)jglob/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13201((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13201[0].decDep();
return;
}
for (int i_darts_counter_temp13201 = (*i);i_darts_counter_temp13201<endRange && i_darts_counter_temp13201<this->inputsTPParent->lastIteration13201;i_darts_counter_temp13201++)
{
{
(*(*iglob)) = (i_darts_counter_temp13201);
{
/*Loop's init*/
(*j) = 0;
int j_darts_counter_temp13201 = (*j);
for(;j_darts_counter_temp13201 < ny;j_darts_counter_temp13201++){
(*(*jglob)) = j_darts_counter_temp13201;
exact((*(*iglob)), (*(*jglob)), 0, &u[(i_darts_counter_temp13201)][j_darts_counter_temp13201][0][0]);
exact((*(*iglob)), (*(*jglob)), nz - 1, &u[(i_darts_counter_temp13201)][j_darts_counter_temp13201][nz - 1][0]);
}
(*j) = j_darts_counter_temp13201;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13201[0].decDep();
}
TP13201::TP13201(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13201** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts13201(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iglob_darts13201(new int*[this->numThreads]),j_darts13201(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jglob_darts13201(new int*[this->numThreads]), initIteration13201(in_initIteration), lastIteration13201(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13201(new _barrierCodelets13201[1]) ,checkInCodelets13202(new _checkInCodelets13202[this->numThreads]){
/*Initialize the loop parameters*/
range13201 = abs (lastIteration13201 - initIteration13201) / 1;
rangePerCodelet13201 = range13201 / numThreads;
minIteration13201 = min<int>(lastIteration13201, initIteration13201);
remainderRange13201 = range13201 % numThreads;
/*Initialize inputs and vars.*/
this->iglob_darts13201 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jglob_darts13201 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13201[0] = _barrierCodelets13201(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13202 * checkInCodelets13202Ptr = (this->checkInCodelets13202);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13202);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13202Ptr) = _checkInCodelets13202(2,1,this,codeletCounter);
#else
(*checkInCodelets13202Ptr) = _checkInCodelets13202(1,1,this,codeletCounter);
#endif
(*checkInCodelets13202Ptr).decDep();
checkInCodelets13202Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13201::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13202[localID].setID (codeletID);
this->checkInCodelets13202[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13202[localID + this->baseNumThreads * i] = _checkInCodelets13202(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13202[localID + this->baseNumThreads * i] = _checkInCodelets13202(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13202[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13202[localID + this->baseNumThreads * i].decDep();
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
TP13201::~TP13201(){
delete [] iglob_darts13201;
delete [] jglob_darts13201;
delete [] barrierCodelets13201;
delete [] checkInCodelets13202;
}
/*TP13252: OMPForDirective*/
void TP13252::_barrierCodelets13252::fire(void)
{
TP13252* myTP = static_cast<TP13252*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13252[0].decDep ();
}
bool TP13252::requestNewRangeIterations13252(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13252 * codeletID;
int tempEndRange   = rangePerCodelet13252 * (codeletID + 1);
if (remainderRange13252 != 0)
{
if (codeletID < (uint32_t)remainderRange13252)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13252;
tempEndRange += remainderRange13252;
}
}
tempStartRange = tempStartRange*1 + minIteration13252;
tempEndRange = tempEndRange*1 + minIteration13252;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13252 < lastIteration13252)
{
(this->inputsTPParent->i_darts13252[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts13252[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13252;
}
}
return isThereNewIteration;
}
void TP13252::_checkInCodelets13253::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iglob_darts13252[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13197[this->getID()]);

/*printing node 13253: ForStmt*/
/*var: i*/
/*var: iglob*/
/*var: k*/
int* i = &(this->inputsTPParent->i_darts13252[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iglob = &(this->inputsTPParent->iglob_darts13252[this->getLocalID()]);
(void)iglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13252[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13252((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13252[0].decDep();
return;
}
for (int i_darts_counter_temp13252 = (*i);i_darts_counter_temp13252<endRange && i_darts_counter_temp13252<this->inputsTPParent->lastIteration13252;i_darts_counter_temp13252++)
{
{
(*(*iglob)) = (i_darts_counter_temp13252);
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp13252 = (*k);
for(;k_darts_counter_temp13252 < nz;k_darts_counter_temp13252++){
exact((*(*iglob)), 0, k_darts_counter_temp13252, &u[(i_darts_counter_temp13252)][0][k_darts_counter_temp13252][0]);
}
(*k) = k_darts_counter_temp13252;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13252[0].decDep();
}
TP13252::TP13252(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13252** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts13252(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iglob_darts13252(new int*[this->numThreads]),k_darts13252(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13252(in_initIteration), lastIteration13252(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13252(new _barrierCodelets13252[1]) ,checkInCodelets13253(new _checkInCodelets13253[this->numThreads]){
/*Initialize the loop parameters*/
range13252 = abs (lastIteration13252 - initIteration13252) / 1;
rangePerCodelet13252 = range13252 / numThreads;
minIteration13252 = min<int>(lastIteration13252, initIteration13252);
remainderRange13252 = range13252 % numThreads;
/*Initialize inputs and vars.*/
this->iglob_darts13252 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13252[0] = _barrierCodelets13252(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13253 * checkInCodelets13253Ptr = (this->checkInCodelets13253);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13253);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13253Ptr) = _checkInCodelets13253(2,1,this,codeletCounter);
#else
(*checkInCodelets13253Ptr) = _checkInCodelets13253(1,1,this,codeletCounter);
#endif
(*checkInCodelets13253Ptr).decDep();
checkInCodelets13253Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13252::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13253[localID].setID (codeletID);
this->checkInCodelets13253[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13253[localID + this->baseNumThreads * i] = _checkInCodelets13253(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13253[localID + this->baseNumThreads * i] = _checkInCodelets13253(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13253[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13253[localID + this->baseNumThreads * i].decDep();
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
TP13252::~TP13252(){
delete [] iglob_darts13252;
delete [] barrierCodelets13252;
delete [] checkInCodelets13253;
}
/*TP13294: OMPForDirective*/
void TP13294::_barrierCodelets13294::fire(void)
{
TP13294* myTP = static_cast<TP13294*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13294[0].decDep ();
}
bool TP13294::requestNewRangeIterations13294(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13294 * codeletID;
int tempEndRange   = rangePerCodelet13294 * (codeletID + 1);
if (remainderRange13294 != 0)
{
if (codeletID < (uint32_t)remainderRange13294)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13294;
tempEndRange += remainderRange13294;
}
}
tempStartRange = tempStartRange*1 + minIteration13294;
tempEndRange = tempEndRange*1 + minIteration13294;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13294 < lastIteration13294)
{
(this->inputsTPParent->i_darts13294[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts13294[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13294;
}
}
return isThereNewIteration;
}
void TP13294::_checkInCodelets13295::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iglob_darts13294[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13197[this->getID()]);

/*printing node 13295: ForStmt*/
/*var: i*/
/*var: iglob*/
/*var: k*/
int* i = &(this->inputsTPParent->i_darts13294[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iglob = &(this->inputsTPParent->iglob_darts13294[this->getLocalID()]);
(void)iglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13294[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13294((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13294[0].decDep();
return;
}
for (int i_darts_counter_temp13294 = (*i);i_darts_counter_temp13294<endRange && i_darts_counter_temp13294<this->inputsTPParent->lastIteration13294;i_darts_counter_temp13294++)
{
{
(*(*iglob)) = (i_darts_counter_temp13294);
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp13294 = (*k);
for(;k_darts_counter_temp13294 < nz;k_darts_counter_temp13294++){
exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13294, &u[(i_darts_counter_temp13294)][ny - 1][k_darts_counter_temp13294][0]);
}
(*k) = k_darts_counter_temp13294;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13294[0].decDep();
}
TP13294::TP13294(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13294** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts13294(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iglob_darts13294(new int*[this->numThreads]),k_darts13294(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13294(in_initIteration), lastIteration13294(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13294(new _barrierCodelets13294[1]) ,checkInCodelets13295(new _checkInCodelets13295[this->numThreads]){
/*Initialize the loop parameters*/
range13294 = abs (lastIteration13294 - initIteration13294) / 1;
rangePerCodelet13294 = range13294 / numThreads;
minIteration13294 = min<int>(lastIteration13294, initIteration13294);
remainderRange13294 = range13294 % numThreads;
/*Initialize inputs and vars.*/
this->iglob_darts13294 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13294[0] = _barrierCodelets13294(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13295 * checkInCodelets13295Ptr = (this->checkInCodelets13295);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13295);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13295Ptr) = _checkInCodelets13295(2,1,this,codeletCounter);
#else
(*checkInCodelets13295Ptr) = _checkInCodelets13295(1,1,this,codeletCounter);
#endif
(*checkInCodelets13295Ptr).decDep();
checkInCodelets13295Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13294::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13295[localID].setID (codeletID);
this->checkInCodelets13295[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13295[localID + this->baseNumThreads * i] = _checkInCodelets13295(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13295[localID + this->baseNumThreads * i] = _checkInCodelets13295(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13295[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13295[localID + this->baseNumThreads * i].decDep();
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
TP13294::~TP13294(){
delete [] iglob_darts13294;
delete [] barrierCodelets13294;
delete [] checkInCodelets13295;
}
/*TP13338: OMPForDirective*/
void TP13338::_barrierCodelets13338::fire(void)
{
TP13338* myTP = static_cast<TP13338*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13338[0].decDep ();
}
bool TP13338::requestNewRangeIterations13338(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13338 * codeletID;
int tempEndRange   = rangePerCodelet13338 * (codeletID + 1);
if (remainderRange13338 != 0)
{
if (codeletID < (uint32_t)remainderRange13338)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13338;
tempEndRange += remainderRange13338;
}
}
tempStartRange = tempStartRange*1 + minIteration13338;
tempEndRange = tempEndRange*1 + minIteration13338;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13338 < lastIteration13338)
{
(this->inputsTPParent->j_darts13338[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->j_darts13338[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13338;
}
}
return isThereNewIteration;
}
void TP13338::_checkInCodelets13339::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jglob_darts13338[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13197[this->getID()]);

/*printing node 13339: ForStmt*/
/*var: j*/
/*var: jglob*/
/*var: k*/
int* j = &(this->inputsTPParent->j_darts13338[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jglob = &(this->inputsTPParent->jglob_darts13338[this->getLocalID()]);
(void)jglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13338[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13338((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13338[0].decDep();
return;
}
for (int j_darts_counter_temp13338 = (*j);j_darts_counter_temp13338<endRange && j_darts_counter_temp13338<this->inputsTPParent->lastIteration13338;j_darts_counter_temp13338++)
{
{
(*(*jglob)) = (j_darts_counter_temp13338);
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp13338 = (*k);
for(;k_darts_counter_temp13338 < nz;k_darts_counter_temp13338++){
exact(0, (*(*jglob)), k_darts_counter_temp13338, &u[0][(j_darts_counter_temp13338)][k_darts_counter_temp13338][0]);
}
(*k) = k_darts_counter_temp13338;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13338[0].decDep();
}
TP13338::TP13338(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13338** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),j_darts13338(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jglob_darts13338(new int*[this->numThreads]),k_darts13338(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13338(in_initIteration), lastIteration13338(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13338(new _barrierCodelets13338[1]) ,checkInCodelets13339(new _checkInCodelets13339[this->numThreads]){
/*Initialize the loop parameters*/
range13338 = abs (lastIteration13338 - initIteration13338) / 1;
rangePerCodelet13338 = range13338 / numThreads;
minIteration13338 = min<int>(lastIteration13338, initIteration13338);
remainderRange13338 = range13338 % numThreads;
/*Initialize inputs and vars.*/
this->jglob_darts13338 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13338[0] = _barrierCodelets13338(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13339 * checkInCodelets13339Ptr = (this->checkInCodelets13339);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13339);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13339Ptr) = _checkInCodelets13339(2,1,this,codeletCounter);
#else
(*checkInCodelets13339Ptr) = _checkInCodelets13339(1,1,this,codeletCounter);
#endif
(*checkInCodelets13339Ptr).decDep();
checkInCodelets13339Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13338::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13339[localID].setID (codeletID);
this->checkInCodelets13339[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13339[localID + this->baseNumThreads * i] = _checkInCodelets13339(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13339[localID + this->baseNumThreads * i] = _checkInCodelets13339(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13339[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13339[localID + this->baseNumThreads * i].decDep();
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
TP13338::~TP13338(){
delete [] jglob_darts13338;
delete [] barrierCodelets13338;
delete [] checkInCodelets13339;
}
/*TP13380: OMPForDirective*/
void TP13380::_barrierCodelets13380::fire(void)
{
TP13380* myTP = static_cast<TP13380*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13380[0].decDep ();
}
bool TP13380::requestNewRangeIterations13380(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13380 * codeletID;
int tempEndRange   = rangePerCodelet13380 * (codeletID + 1);
if (remainderRange13380 != 0)
{
if (codeletID < (uint32_t)remainderRange13380)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13380;
tempEndRange += remainderRange13380;
}
}
tempStartRange = tempStartRange*1 + minIteration13380;
tempEndRange = tempEndRange*1 + minIteration13380;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13380 < lastIteration13380)
{
(this->inputsTPParent->j_darts13380[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->j_darts13380[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13380;
}
}
return isThereNewIteration;
}
void TP13380::_checkInCodelets13381::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jglob_darts13380[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13197[this->getID()]);

/*printing node 13381: ForStmt*/
/*var: j*/
/*var: jglob*/
/*var: k*/
int* j = &(this->inputsTPParent->j_darts13380[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jglob = &(this->inputsTPParent->jglob_darts13380[this->getLocalID()]);
(void)jglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13380[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13380((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13380[0].decDep();
return;
}
for (int j_darts_counter_temp13380 = (*j);j_darts_counter_temp13380<endRange && j_darts_counter_temp13380<this->inputsTPParent->lastIteration13380;j_darts_counter_temp13380++)
{
{
(*(*jglob)) = (j_darts_counter_temp13380);
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp13380 = (*k);
for(;k_darts_counter_temp13380 < nz;k_darts_counter_temp13380++){
exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13380, &u[nx - 1][(j_darts_counter_temp13380)][k_darts_counter_temp13380][0]);
}
(*k) = k_darts_counter_temp13380;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13380[0].decDep();
}
TP13380::TP13380(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13380** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),j_darts13380(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jglob_darts13380(new int*[this->numThreads]),k_darts13380(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13380(in_initIteration), lastIteration13380(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13380(new _barrierCodelets13380[1]) ,checkInCodelets13381(new _checkInCodelets13381[this->numThreads]){
/*Initialize the loop parameters*/
range13380 = abs (lastIteration13380 - initIteration13380) / 1;
rangePerCodelet13380 = range13380 / numThreads;
minIteration13380 = min<int>(lastIteration13380, initIteration13380);
remainderRange13380 = range13380 % numThreads;
/*Initialize inputs and vars.*/
this->jglob_darts13380 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13380[0] = _barrierCodelets13380(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13381 * checkInCodelets13381Ptr = (this->checkInCodelets13381);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13381);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13381Ptr) = _checkInCodelets13381(2,1,this,codeletCounter);
#else
(*checkInCodelets13381Ptr) = _checkInCodelets13381(1,1,this,codeletCounter);
#endif
(*checkInCodelets13381Ptr).decDep();
checkInCodelets13381Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13380::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13381[localID].setID (codeletID);
this->checkInCodelets13381[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13381[localID + this->baseNumThreads * i] = _checkInCodelets13381(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13381[localID + this->baseNumThreads * i] = _checkInCodelets13381(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13381[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13381[localID + this->baseNumThreads * i].decDep();
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
TP13380::~TP13380(){
delete [] jglob_darts13380;
delete [] barrierCodelets13380;
delete [] checkInCodelets13381;
}
/*TP13764: OMPParallelDirective*/
void TP13764::_barrierCodelets13764::fire(void)
{
TP13764* myTP = static_cast<TP13764*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP13764::_checkInCodelets13770::fire(void)
{
/*region 13770 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13770;
if(idx < myTP->TPsToUse13770){
if (!__sync_val_compare_and_swap (&(myTP->TP13770_alreadyLaunched[idx]), 0, 1)){
int range = abs (ny - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13770;
int minIteration = min<int >(ny, 0);
int remainderRange = range % myTP->TPsToUse13770;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < ny)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13770 - 1)
{
lastIteration = ny;
}
#if USEINVOKE == 1
invoke < TP13770 > (myTP, myTP->codeletsPerTP13770 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13770Ptr[idx]));
#else
place < TP13770 > (idx, myTP, myTP->codeletsPerTP13770 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13770Ptr[idx]));
#endif
}else{
if (myTP->TP13770Ptr[idx] != nullptr){
myTP->TP13770Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13764::_barrierCodelets13770::fire(void)
{
TP13764* myTP =  static_cast<TP13764*>(myTP_);
myTP->TPParent->barrierCodelets13764[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets13764[0]));
}
TP13764::TP13764(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),eta_darts13764(new double[this->numThreads])/*VARIABLE*/,i_darts13764(new int[this->numThreads])/*VARIABLE*/,iglob_darts13764(new int[this->numThreads])/*VARIABLE*/,j_darts13764(new int[this->numThreads])/*VARIABLE*/,jglob_darts13764(new int[this->numThreads])/*VARIABLE*/,k_darts13764(new int[this->numThreads])/*VARIABLE*/,m_darts13764(new int[this->numThreads])/*VARIABLE*/,peta_darts13764(new double[this->numThreads])/*VARIABLE*/,pxi_darts13764(new double[this->numThreads])/*VARIABLE*/,pzeta_darts13764(new double[this->numThreads])/*VARIABLE*/,xi_darts13764(new double[this->numThreads])/*VARIABLE*/,zeta_darts13764(new double[this->numThreads])/*VARIABLE*/, TP13770Ptr(new TP13770 *[NUMTPS13770]), TP13770_alreadyLaunched(new size_t [NUMTPS13770]), numTPsSet13770(0), numTPsReady13770(0), TPsToUse13770(NUMTPS13770), codeletsPerTP13770(this->numThreads/NUMTPS13770), totalCodelets13770(this->TPsToUse13770*this->codeletsPerTP13770) ,barrierCodelets13764(new _barrierCodelets13764[1]) ,checkInCodelets13770(new _checkInCodelets13770[this->numThreads]) ,barrierCodelets13770(new _barrierCodelets13770[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13764[0] = _barrierCodelets13764(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets13770[0] = _barrierCodelets13770(NUMTPS13770,NUMTPS13770,this, 0);
_checkInCodelets13770 * checkInCodelets13770Ptr = (this->checkInCodelets13770);
for(int i=0; i<NUMTPS13770; i++)
{
TP13770Ptr[i] = nullptr;
TP13770_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets13770Ptr) = _checkInCodelets13770(1,1,this,codeletCounter);
(*checkInCodelets13770Ptr).decDep();
checkInCodelets13770Ptr++;
}
}
TP13764::~TP13764(){
delete []eta_darts13764;
delete []i_darts13764;
delete []iglob_darts13764;
delete []j_darts13764;
delete []jglob_darts13764;
delete []k_darts13764;
delete []m_darts13764;
delete []peta_darts13764;
delete []pxi_darts13764;
delete []pzeta_darts13764;
delete []xi_darts13764;
delete []zeta_darts13764;
delete [] barrierCodelets13764;
delete [] barrierCodelets13770;
delete [] checkInCodelets13770;
}
/*TP13770: OMPForDirective*/
void TP13770::_barrierCodelets13770::fire(void)
{
TP13770* myTP = static_cast<TP13770*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13770[0].decDep ();
}
bool TP13770::requestNewRangeIterations13770(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13770 * codeletID;
int tempEndRange   = rangePerCodelet13770 * (codeletID + 1);
if (remainderRange13770 != 0)
{
if (codeletID < (uint32_t)remainderRange13770)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13770;
tempEndRange += remainderRange13770;
}
}
tempStartRange = tempStartRange*1 + minIteration13770;
tempEndRange = tempEndRange*1 + minIteration13770;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13770 < lastIteration13770)
{
(this->inputsTPParent->j_darts13770[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->j_darts13770[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13770;
}
}
return isThereNewIteration;
}
void TP13770::_checkInCodelets13771::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->eta_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->eta_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->iglob_darts13770[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->jglob_darts13770[this->getLocalID()] = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->peta_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->peta_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->pxi_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->pxi_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->pzeta_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->pzeta_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->xi_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->xi_darts13764[this->getID()]);
/*Get pointer from parent for variable
 with shared scope in this region but private
 in the enclosing one.*/
this->inputsTPParent->zeta_darts13770[this->getLocalID()] = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts13764[this->getID()]);

/*printing node 13771: ForStmt*/
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
double** eta = &(this->inputsTPParent->eta_darts13770[this->getLocalID()]);
(void)eta/*OMP_SHARED_PRIVATE*/;
int* i = &(this->inputsTPParent->i_darts13770[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int** iglob = &(this->inputsTPParent->iglob_darts13770[this->getLocalID()]);
(void)iglob/*OMP_SHARED_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts13770[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int** jglob = &(this->inputsTPParent->jglob_darts13770[this->getLocalID()]);
(void)jglob/*OMP_SHARED_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13770[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts13770[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double** peta = &(this->inputsTPParent->peta_darts13770[this->getLocalID()]);
(void)peta/*OMP_SHARED_PRIVATE*/;
double** pxi = &(this->inputsTPParent->pxi_darts13770[this->getLocalID()]);
(void)pxi/*OMP_SHARED_PRIVATE*/;
double** pzeta = &(this->inputsTPParent->pzeta_darts13770[this->getLocalID()]);
(void)pzeta/*OMP_SHARED_PRIVATE*/;
double** xi = &(this->inputsTPParent->xi_darts13770[this->getLocalID()]);
(void)xi/*OMP_SHARED_PRIVATE*/;
double** zeta = &(this->inputsTPParent->zeta_darts13770[this->getLocalID()]);
(void)zeta/*OMP_SHARED_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13770((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13770[0].decDep();
return;
}
for (int j_darts_counter_temp13770 = (*j);j_darts_counter_temp13770<endRange && j_darts_counter_temp13770<this->inputsTPParent->lastIteration13770;j_darts_counter_temp13770++)
{
{
(*(*jglob)) = (j_darts_counter_temp13770);
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp13770 = (*k);
for(;k_darts_counter_temp13770 < nz - 1;k_darts_counter_temp13770++){
(*(*zeta)) = ((double)k_darts_counter_temp13770) / (nz - 1);
if((*(*jglob)) != 0 && (*(*jglob)) != ny0 - 1)
{
(*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
{
/*Loop's init*/
(*i) = 0;
int i_darts_counter_temp13770 = (*i);
for(;i_darts_counter_temp13770 < nx;i_darts_counter_temp13770++){
(*(*iglob)) = i_darts_counter_temp13770;
if((*(*iglob)) != 0 && (*(*iglob)) != nx0 - 1)
{
(*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
exact(0, (*(*jglob)), k_darts_counter_temp13770, ue_1jk);
exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13770, ue_nx0jk);
exact((*(*iglob)), 0, k_darts_counter_temp13770, ue_i1k);
exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13770, ue_iny0k);
exact((*(*iglob)), (*(*jglob)), 0, ue_ij1);
exact((*(*iglob)), (*(*jglob)), nz - 1, ue_ijnz);
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp13770 = (*m);
for(;m_darts_counter_temp13770 < 5;m_darts_counter_temp13770++){
(*(*pxi)) = (1. - (*(*xi))) * ue_1jk[m_darts_counter_temp13770] + (*(*xi)) * ue_nx0jk[m_darts_counter_temp13770];
(*(*peta)) = (1. - (*(*eta))) * ue_i1k[m_darts_counter_temp13770] + (*(*eta)) * ue_iny0k[m_darts_counter_temp13770];
(*(*pzeta)) = (1. - (*(*zeta))) * ue_ij1[m_darts_counter_temp13770] + (*(*zeta)) * ue_ijnz[m_darts_counter_temp13770];
u[i_darts_counter_temp13770][(j_darts_counter_temp13770)][k_darts_counter_temp13770][m_darts_counter_temp13770] = (*(*pxi)) + (*(*peta)) + (*(*pzeta)) - (*(*pxi)) * (*(*peta)) - (*(*peta)) * (*(*pzeta)) - (*(*pzeta)) * (*(*pxi)) + (*(*pxi)) * (*(*peta)) * (*(*pzeta));
}
(*m) = m_darts_counter_temp13770;
}
}
}
(*i) = i_darts_counter_temp13770;
}
}
}
(*k) = k_darts_counter_temp13770;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13770[0].decDep();
}
TP13770::TP13770(int in_numThreads, int in_mainCodeletID, TP13764* in_TPParent, int in_initIteration, int in_lastIteration, TP13770** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),eta_darts13770(new double*[this->numThreads]),i_darts13770(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,iglob_darts13770(new int*[this->numThreads]),j_darts13770(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,jglob_darts13770(new int*[this->numThreads]),k_darts13770(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts13770(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,peta_darts13770(new double*[this->numThreads]),pxi_darts13770(new double*[this->numThreads]),pzeta_darts13770(new double*[this->numThreads]),xi_darts13770(new double*[this->numThreads]),zeta_darts13770(new double*[this->numThreads]), initIteration13770(in_initIteration), lastIteration13770(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13770(new _barrierCodelets13770[1]) ,checkInCodelets13771(new _checkInCodelets13771[this->numThreads]){
/*Initialize the loop parameters*/
range13770 = abs (lastIteration13770 - initIteration13770) / 1;
rangePerCodelet13770 = range13770 / numThreads;
minIteration13770 = min<int>(lastIteration13770, initIteration13770);
remainderRange13770 = range13770 % numThreads;
/*Initialize inputs and vars.*/
this->eta_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->iglob_darts13770 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->jglob_darts13770 = (int**)malloc(sizeof(int*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->peta_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->pxi_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->pzeta_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->xi_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
this->zeta_darts13770 = (double**)malloc(sizeof(double*) * this->numThreads)/*OMP_SHARED_PRIVATE*/;
/*Initialize Codelets*/
barrierCodelets13770[0] = _barrierCodelets13770(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13771 * checkInCodelets13771Ptr = (this->checkInCodelets13771);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13771);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13771Ptr) = _checkInCodelets13771(2,1,this,codeletCounter);
#else
(*checkInCodelets13771Ptr) = _checkInCodelets13771(1,1,this,codeletCounter);
#endif
(*checkInCodelets13771Ptr).decDep();
checkInCodelets13771Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13770::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13771[localID].setID (codeletID);
this->checkInCodelets13771[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13771[localID + this->baseNumThreads * i] = _checkInCodelets13771(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13771[localID + this->baseNumThreads * i] = _checkInCodelets13771(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13771[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13771[localID + this->baseNumThreads * i].decDep();
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
TP13770::~TP13770(){
delete [] eta_darts13770;
delete [] iglob_darts13770;
delete [] jglob_darts13770;
delete [] peta_darts13770;
delete [] pxi_darts13770;
delete [] pzeta_darts13770;
delete [] xi_darts13770;
delete [] zeta_darts13770;
delete [] barrierCodelets13770;
delete [] checkInCodelets13771;
}
/*TP13896: OMPParallelDirective*/
void TP13896::_barrierCodelets13896::fire(void)
{
TP13896* myTP = static_cast<TP13896*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP13896::_checkInCodelets13898::fire(void)
{
/*region 13898 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13898;
if(idx < myTP->TPsToUse13898){
if (!__sync_val_compare_and_swap (&(myTP->TP13898_alreadyLaunched[idx]), 0, 1)){
int range = abs (12 - 0) / 1;
int rangePerCodelet = range / myTP->TPsToUse13898;
int minIteration = min<int >(12, 0);
int remainderRange = range % myTP->TPsToUse13898;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(0 < 12)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13898 - 1)
{
lastIteration = 12;
}
#if USEINVOKE == 1
invoke < TP13898 > (myTP, myTP->codeletsPerTP13898 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13898Ptr[idx]));
#else
place < TP13898 > (idx, myTP, myTP->codeletsPerTP13898 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13898Ptr[idx]));
#endif
}else{
if (myTP->TP13898Ptr[idx] != nullptr){
myTP->TP13898Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13896::_barrierCodelets13898::fire(void)
{
TP13896* myTP =  static_cast<TP13896*>(myTP_);
myTP->TPParent->barrierCodelets13896[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets13896[0]));
}
TP13896::TP13896(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),i_darts13896(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts13896(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts13896(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts13896(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, TP13898Ptr(new TP13898 *[NUMTPS13898]), TP13898_alreadyLaunched(new size_t [NUMTPS13898]), numTPsSet13898(0), numTPsReady13898(0), TPsToUse13898(NUMTPS13898), codeletsPerTP13898(this->numThreads/NUMTPS13898), totalCodelets13898(this->TPsToUse13898*this->codeletsPerTP13898) ,barrierCodelets13896(new _barrierCodelets13896[1]) ,checkInCodelets13898(new _checkInCodelets13898[this->numThreads]) ,barrierCodelets13898(new _barrierCodelets13898[1]){
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13896[0] = _barrierCodelets13896(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets13898[0] = _barrierCodelets13898(NUMTPS13898,NUMTPS13898,this, 0);
_checkInCodelets13898 * checkInCodelets13898Ptr = (this->checkInCodelets13898);
for(int i=0; i<NUMTPS13898; i++)
{
TP13898Ptr[i] = nullptr;
TP13898_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets13898Ptr) = _checkInCodelets13898(1,1,this,codeletCounter);
(*checkInCodelets13898Ptr).decDep();
checkInCodelets13898Ptr++;
}
}
TP13896::~TP13896(){
delete [] barrierCodelets13896;
delete [] barrierCodelets13898;
delete [] checkInCodelets13898;
}
/*TP13898: OMPForDirective*/
void TP13898::_barrierCodelets13898::fire(void)
{
TP13898* myTP = static_cast<TP13898*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13898[0].decDep ();
}
bool TP13898::requestNewRangeIterations13898(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13898 * codeletID;
int tempEndRange   = rangePerCodelet13898 * (codeletID + 1);
if (remainderRange13898 != 0)
{
if (codeletID < (uint32_t)remainderRange13898)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13898;
tempEndRange += remainderRange13898;
}
}
tempStartRange = tempStartRange*1 + minIteration13898;
tempEndRange = tempEndRange*1 + minIteration13898;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13898 < lastIteration13898)
{
(this->inputsTPParent->i_darts13898[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts13898[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13898;
}
}
return isThereNewIteration;
}
void TP13898::_checkInCodelets13899::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 13899: ForStmt*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: m*/
int* i = &(this->inputsTPParent->i_darts13898[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts13898[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13898[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts13898[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13898((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13898[0].decDep();
return;
}
for (int i_darts_counter_temp13898 = (*i);i_darts_counter_temp13898<endRange && i_darts_counter_temp13898<this->inputsTPParent->lastIteration13898;i_darts_counter_temp13898++)
{
{
{
/*Loop's init*/
(*j) = 0;
int j_darts_counter_temp13898 = (*j);
for(;j_darts_counter_temp13898 < 12;j_darts_counter_temp13898++){
{
/*Loop's init*/
(*k) = 0;
int k_darts_counter_temp13898 = (*k);
for(;k_darts_counter_temp13898 < 5;k_darts_counter_temp13898++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp13898 = (*m);
for(;m_darts_counter_temp13898 < 5;m_darts_counter_temp13898++){
a[(i_darts_counter_temp13898)][j_darts_counter_temp13898][k_darts_counter_temp13898][m_darts_counter_temp13898] = 0.;
b[(i_darts_counter_temp13898)][j_darts_counter_temp13898][k_darts_counter_temp13898][m_darts_counter_temp13898] = 0.;
c[(i_darts_counter_temp13898)][j_darts_counter_temp13898][k_darts_counter_temp13898][m_darts_counter_temp13898] = 0.;
d[(i_darts_counter_temp13898)][j_darts_counter_temp13898][k_darts_counter_temp13898][m_darts_counter_temp13898] = 0.;
}
(*m) = m_darts_counter_temp13898;
}
}
(*k) = k_darts_counter_temp13898;
}
}
(*j) = j_darts_counter_temp13898;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13898[0].decDep();
}
TP13898::TP13898(int in_numThreads, int in_mainCodeletID, TP13896* in_TPParent, int in_initIteration, int in_lastIteration, TP13898** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts13898(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts13898(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts13898(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts13898(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13898(in_initIteration), lastIteration13898(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13898(new _barrierCodelets13898[1]) ,checkInCodelets13899(new _checkInCodelets13899[this->numThreads]){
/*Initialize the loop parameters*/
range13898 = abs (lastIteration13898 - initIteration13898) / 1;
rangePerCodelet13898 = range13898 / numThreads;
minIteration13898 = min<int>(lastIteration13898, initIteration13898);
remainderRange13898 = range13898 % numThreads;
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13898[0] = _barrierCodelets13898(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13899 * checkInCodelets13899Ptr = (this->checkInCodelets13899);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13899);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13899Ptr) = _checkInCodelets13899(2,1,this,codeletCounter);
#else
(*checkInCodelets13899Ptr) = _checkInCodelets13899(1,1,this,codeletCounter);
#endif
(*checkInCodelets13899Ptr).decDep();
checkInCodelets13899Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13898::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13899[localID].setID (codeletID);
this->checkInCodelets13899[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13899[localID + this->baseNumThreads * i] = _checkInCodelets13899(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13899[localID + this->baseNumThreads * i] = _checkInCodelets13899(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13899[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13899[localID + this->baseNumThreads * i].decDep();
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
TP13898::~TP13898(){
delete [] barrierCodelets13898;
delete [] checkInCodelets13899;
}
/*TP13995: OMPParallelDirective*/
void TP13995::_barrierCodelets13995::fire(void)
{
TP13995* myTP = static_cast<TP13995*>(myTP_);
myTP->controlTPParent->nextCodelet->decDep ();
}
void TP13995::_checkInCodelets13997::fire(void)
{
/*region 13997 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP13997;
if(idx < myTP->TPsToUse13997){
if (!__sync_val_compare_and_swap (&(myTP->TP13997_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse13997;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse13997;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse13997 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse13997 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP13997 > (myTP, myTP->codeletsPerTP13997 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13997Ptr[idx]));
#else
place < TP13997 > (idx, myTP, myTP->codeletsPerTP13997 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(myTP->TP13997Ptr[idx]));
#endif
}else{
if (myTP->TP13997Ptr[idx] != nullptr){
myTP->TP13997Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13995::_barrierCodelets13997::fire(void)
{
TP13995* myTP =  static_cast<TP13995*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets14052[codeletsCounter].decDep();
}
}
}
void TP13995::_checkInCodelets14052::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/

/*printing node 14053: CallExpr*/
printf("init complete...\n");
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 14052 nextRegion: 14055 */
myTP->controlTPParent->checkInCodelets14055[this->getID()].decDep();
}
else
{
/*Signaling next codelet region: 14052 nextRegion: 14055 */
myTP->checkInCodelets14055[this->getID()].decDep();
}
}
void TP13995::_checkInCodelets14055::fire(void)
{

/*printing node 14055: BinaryOperator*/
(this->inputsTPParent->k_darts13995[this->getID()]) = 1;

/*printing node 14056: BinaryOperator*/
/*Print the code for a condition node in a complex loop stmt */
if((this->inputsTPParent->k_darts13995[this->getID()]) <= nz - 2){
/*Signal the first codelet in the loop*/
myTP->checkInCodelets14054[this->getID()].decDep();
return;
}
else{
/*Signal the codelet after the loop from the end condional node.*/
/*Signaling next codelet region: 14058 nextRegion: 14066 */
myTP->controlTPParent->barrierCodelets14066[0].decDep();
return;
}
}
void TP13995::_checkInCodelets14054::fire(void)
{

/*printing node 14054: ForStmt*/
bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14054_LoopCounter), myTP->controlTPParent->TP14054_LoopCounterPerThread[this->getID ()], myTP->controlTPParent->TP14054_LoopCounterPerThread[this->getID ()] + 1);
unsigned int iterIdx = myTP->controlTPParent->TP14054_LoopCounterPerThread[this->getID ()];
if (haveToLaunch)
{
this->resetCodelet(); 
myTP->controlTPParent->TP14054PtrVec.push_back(nullptr);
myTP->controlTPParent->TP14054_LoopCounterPerThread[this->getID ()] += 1;
invoke < TP14054 > (myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent, &(myTP->controlTPParent->TP14054PtrVec.back()));
}
else
{
if(myTP->controlTPParent->TP14054PtrVec.size () == 0)
{
this->resetCodelet ();
this->decDep ();
return;
} else if (myTP->controlTPParent->TP14054PtrVec.size () < (iterIdx + 1))
{
this->resetCodelet ();
this->decDep ();
return;
} else if (myTP->controlTPParent->TP14054PtrVec[iterIdx] == nullptr)
{
this->resetCodelet ();
this->decDep ();
return;
}
else
{
this->resetCodelet ();
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14054PtrVec[iterIdx]->firstCodelet[this->getID ()].decDep();
#else
myTP->controlTPParent->TP14054PtrVec[iterIdx]->availableCodelets[this->getID ()] = 1;
#endif
myTP->controlTPParent->TP14054_LoopCounterPerThread[this->getID ()] += 1;
}
}
}
void TP13995::_checkInCodelets14058::fire(void)
{

/*printing node 14058: UnaryOperator*/
(this->inputsTPParent->k_darts13995[this->getID()])++;

/*printing node 14534: BinaryOperator*/
/*Print the code for a condition node in a complex loop stmt */
if((this->inputsTPParent->k_darts13995[this->getID()]) <= nz - 2){
this->resetCodelet();
/*Signal the first codelet in the loop*/
myTP->checkInCodelets14054[this->getID()].decDep();
return;
}
else{
/*Signal the codelet after the loop from the condtional node.*/
/*Signaling next codelet region: 14058 nextRegion: 14066 */
myTP->controlTPParent->barrierCodelets14066[0].decDep();
return;
}
}
void TP13995::_barrierCodelets14066::fire(void)
{
TP13995* myTP =  static_cast<TP13995*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets14067[codeletsCounter].decDep();
}
}
}
void TP13995::_checkInCodelets14067::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/

/*printing node 14068: CallExpr*/
printf("lower triangular part complete...\n");
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 14067 nextRegion: 14070 */
myTP->controlTPParent->checkInCodelets14070[this->getID()].decDep();
}
else
{
/*Signaling next codelet region: 14067 nextRegion: 14070 */
myTP->checkInCodelets14070[this->getID()].decDep();
}
}
void TP13995::_checkInCodelets14070::fire(void)
{

/*printing node 14070: BinaryOperator*/
(this->inputsTPParent->k_darts13995[this->getID()]) = nz - 2;

/*printing node 14072: BinaryOperator*/
/*Print the code for a condition node in a complex loop stmt */
if((this->inputsTPParent->k_darts13995[this->getID()]) >= 1){
/*Signal the first codelet in the loop*/
myTP->checkInCodelets14069[this->getID()].decDep();
return;
}
else{
/*Signal the codelet after the loop from the end condional node.*/
/*Signaling next codelet region: 14073 nextRegion: 14077 */
myTP->controlTPParent->barrierCodelets14077[0].decDep();
return;
}
}
void TP13995::_checkInCodelets14069::fire(void)
{

/*printing node 14069: ForStmt*/
bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14069_LoopCounter), myTP->controlTPParent->TP14069_LoopCounterPerThread[this->getID ()], myTP->controlTPParent->TP14069_LoopCounterPerThread[this->getID ()] + 1);
unsigned int iterIdx = myTP->controlTPParent->TP14069_LoopCounterPerThread[this->getID ()];
if (haveToLaunch)
{
this->resetCodelet(); 
myTP->controlTPParent->TP14069PtrVec.push_back(nullptr);
myTP->controlTPParent->TP14069_LoopCounterPerThread[this->getID ()] += 1;
invoke < TP14069 > (myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent, &(myTP->controlTPParent->TP14069PtrVec.back()));
}
else
{
if(myTP->controlTPParent->TP14069PtrVec.size () == 0)
{
this->resetCodelet ();
this->decDep ();
return;
} else if (myTP->controlTPParent->TP14069PtrVec.size () < (iterIdx + 1))
{
this->resetCodelet ();
this->decDep ();
return;
} else if (myTP->controlTPParent->TP14069PtrVec[iterIdx] == nullptr)
{
this->resetCodelet ();
this->decDep ();
return;
}
else
{
this->resetCodelet ();
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14069PtrVec[iterIdx]->firstCodelet[this->getID ()].decDep();
#else
myTP->controlTPParent->TP14069PtrVec[iterIdx]->availableCodelets[this->getID ()] = 1;
#endif
myTP->controlTPParent->TP14069_LoopCounterPerThread[this->getID ()] += 1;
}
}
}
void TP13995::_checkInCodelets14073::fire(void)
{

/*printing node 14073: UnaryOperator*/
(this->inputsTPParent->k_darts13995[this->getID()])--;

/*printing node 14535: BinaryOperator*/
/*Print the code for a condition node in a complex loop stmt */
if((this->inputsTPParent->k_darts13995[this->getID()]) >= 1){
this->resetCodelet();
/*Signal the first codelet in the loop*/
myTP->checkInCodelets14069[this->getID()].decDep();
return;
}
else{
/*Signal the codelet after the loop from the condtional node.*/
/*Signaling next codelet region: 14073 nextRegion: 14077 */
myTP->controlTPParent->barrierCodelets14077[0].decDep();
return;
}
}
void TP13995::_barrierCodelets14077::fire(void)
{
TP13995* myTP =  static_cast<TP13995*>(myTP_);
{
for(size_t codeletsCounter=0; codeletsCounter < myTP->numThreads;codeletsCounter++)
{
myTP->checkInCodelets14078[codeletsCounter].decDep();
}
}
}
void TP13995::_checkInCodelets14078::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/

/*printing node 14079: CallExpr*/
printf("buts complete...\n");
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 14078 nextRegion: 14080 */
myTP->controlTPParent->checkInCodelets14080[this->getID()].decDep();
}
else
{
/*Signaling next codelet region: 14078 nextRegion: 14080 */
myTP->checkInCodelets14080[this->getID()].decDep();
}
}
void TP13995::_checkInCodelets14080::fire(void)
{
/*region 14080 0*/
/*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if it signals it using dispatchCodelet()*/
size_t idx = this->getID() / myTP->codeletsPerTP14080;
if(idx < myTP->TPsToUse14080){
if (!__sync_val_compare_and_swap (&(myTP->TP14080_alreadyLaunched[idx]), 0, 1)){
int range = abs (iend - ist) / 1;
int rangePerCodelet = range / myTP->TPsToUse14080;
int minIteration = min<int >(iend, ist);
int remainderRange = range % myTP->TPsToUse14080;
int initIteration = rangePerCodelet * idx;
int lastIteration = rangePerCodelet * (idx + 1);
if (remainderRange != 0)
{
if (idx < (uint32_t)remainderRange)
{
initIteration += idx;
lastIteration += (idx + 1);
}
else
{
initIteration += remainderRange;
lastIteration += remainderRange;
}
}
initIteration = initIteration*1 + minIteration;
lastIteration = lastIteration*1 + minIteration;
if(ist < iend)
{
initIteration = min(initIteration, lastIteration);
lastIteration = max(initIteration, lastIteration);
}
else
{
initIteration = max(initIteration, lastIteration);
lastIteration = min(initIteration, lastIteration);
}
if(idx == myTP->TPsToUse14080 - 1)
{
lastIteration = lastIteration + 1;
}
if(idx == myTP->TPsToUse14080 - 1)
{
lastIteration = iend;
}
#if USEINVOKE == 1
invoke < TP14080 > (myTP, myTP->codeletsPerTP14080 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts13995)), &(myTP->TP14080Ptr[idx]));
#else
place < TP14080 > (idx, myTP, myTP->codeletsPerTP14080 * DARTS_CODELETS_MULT , this->getID(), myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts13995)), &(myTP->TP14080Ptr[idx]));
#endif
}else{
if (myTP->TP14080Ptr[idx] != nullptr){
myTP->TP14080Ptr[idx]->dispatchCodelet(this->getID ());
}else{
this->resetCodelet ();
this->decDep ();
}
}
}
}
void TP13995::_barrierCodelets14080::fire(void)
{
TP13995* myTP =  static_cast<TP13995*>(myTP_);
myTP->TPParent->barrierCodelets13995[0].setDep(0);
myTP->add(&(myTP->TPParent->barrierCodelets13995[0]));
}
TP13995::TP13995(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, double *in_tmp):ThreadedProcedure(in_numThreads, in_mainCodeletID), nextCodelet(in_nextCodelet), TPParent(this), controlTPParent(this), inputsTPParent(this),i_darts13995(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,istep_darts13995(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts13995(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts13995(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts13995(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts13995(in_tmp)/*OMP_SHARED - INPUT*/, TP13997Ptr(new TP13997 *[NUMTPS13997]), TP13997_alreadyLaunched(new size_t [NUMTPS13997]), numTPsSet13997(0), numTPsReady13997(0), TPsToUse13997(NUMTPS13997), codeletsPerTP13997(this->numThreads/NUMTPS13997), totalCodelets13997(this->TPsToUse13997*this->codeletsPerTP13997), TP14052_alreadyLaunched(0), TP14054_LoopCounter(0), TP14054_LoopCounterPerThread(new unsigned int[this->numThreads]), TP14067_alreadyLaunched(0), TP14069_LoopCounter(0), TP14069_LoopCounterPerThread(new unsigned int[this->numThreads]), TP14078_alreadyLaunched(0), TP14080Ptr(new TP14080 *[NUMTPS14080]), TP14080_alreadyLaunched(new size_t [NUMTPS14080]), numTPsSet14080(0), numTPsReady14080(0), TPsToUse14080(NUMTPS14080), codeletsPerTP14080(this->numThreads/NUMTPS14080), totalCodelets14080(this->TPsToUse14080*this->codeletsPerTP14080) ,barrierCodelets13995(new _barrierCodelets13995[1]) ,checkInCodelets13997(new _checkInCodelets13997[this->numThreads]) ,barrierCodelets13997(new _barrierCodelets13997[1]) ,checkInCodelets14052(new _checkInCodelets14052[this->numThreads]) ,checkInCodelets14055(new _checkInCodelets14055[this->numThreads]) ,checkInCodelets14054(new _checkInCodelets14054[this->numThreads]) ,checkInCodelets14058(new _checkInCodelets14058[this->numThreads]) ,barrierCodelets14066(new _barrierCodelets14066[1]) ,checkInCodelets14067(new _checkInCodelets14067[this->numThreads]) ,checkInCodelets14070(new _checkInCodelets14070[this->numThreads]) ,checkInCodelets14069(new _checkInCodelets14069[this->numThreads]) ,checkInCodelets14073(new _checkInCodelets14073[this->numThreads]) ,barrierCodelets14077(new _barrierCodelets14077[1]) ,checkInCodelets14078(new _checkInCodelets14078[this->numThreads]) ,checkInCodelets14080(new _checkInCodelets14080[this->numThreads]) ,barrierCodelets14080(new _barrierCodelets14080[1]){
memset((void*)TP14054_LoopCounterPerThread, 0, this->numThreads*sizeof(unsigned int));
memset((void*)TP14069_LoopCounterPerThread, 0, this->numThreads*sizeof(unsigned int));
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13995[0] = _barrierCodelets13995(ompNumThreads,ompNumThreads,this, 0);
barrierCodelets14080[0] = _barrierCodelets14080(NUMTPS14080,NUMTPS14080,this, 0);
barrierCodelets14077[0] = _barrierCodelets14077(this->numThreads,this->numThreads,this, 0);
barrierCodelets14066[0] = _barrierCodelets14066(this->numThreads,this->numThreads,this, 0);
barrierCodelets13997[0] = _barrierCodelets13997(NUMTPS13997,NUMTPS13997,this, 0);
_checkInCodelets14080 * checkInCodelets14080Ptr = (this->checkInCodelets14080);
for(int i=0; i<NUMTPS14080; i++)
{
TP14080Ptr[i] = nullptr;
TP14080_alreadyLaunched[i] = 0;
}
_checkInCodelets14078 * checkInCodelets14078Ptr = (this->checkInCodelets14078);
_checkInCodelets14073 * checkInCodelets14073Ptr = (this->checkInCodelets14073);
_checkInCodelets14069 * checkInCodelets14069Ptr = (this->checkInCodelets14069);
_checkInCodelets14070 * checkInCodelets14070Ptr = (this->checkInCodelets14070);
_checkInCodelets14067 * checkInCodelets14067Ptr = (this->checkInCodelets14067);
_checkInCodelets14058 * checkInCodelets14058Ptr = (this->checkInCodelets14058);
_checkInCodelets14054 * checkInCodelets14054Ptr = (this->checkInCodelets14054);
_checkInCodelets14055 * checkInCodelets14055Ptr = (this->checkInCodelets14055);
_checkInCodelets14052 * checkInCodelets14052Ptr = (this->checkInCodelets14052);
_checkInCodelets13997 * checkInCodelets13997Ptr = (this->checkInCodelets13997);
for(int i=0; i<NUMTPS13997; i++)
{
TP13997Ptr[i] = nullptr;
TP13997_alreadyLaunched[i] = 0;
}
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets14080Ptr) = _checkInCodelets14080(1,1,this,codeletCounter);
checkInCodelets14080Ptr++;
(*checkInCodelets14078Ptr) = _checkInCodelets14078(1,1,this,codeletCounter);
checkInCodelets14078Ptr++;
(*checkInCodelets14073Ptr) = _checkInCodelets14073(1,1,this,codeletCounter);
checkInCodelets14073Ptr++;
(*checkInCodelets14069Ptr) = _checkInCodelets14069(1,1,this,codeletCounter);
checkInCodelets14069Ptr++;
(*checkInCodelets14070Ptr) = _checkInCodelets14070(1,1,this,codeletCounter);
checkInCodelets14070Ptr++;
(*checkInCodelets14067Ptr) = _checkInCodelets14067(1,1,this,codeletCounter);
checkInCodelets14067Ptr++;
(*checkInCodelets14058Ptr) = _checkInCodelets14058(1,1,this,codeletCounter);
checkInCodelets14058Ptr++;
(*checkInCodelets14054Ptr) = _checkInCodelets14054(1,1,this,codeletCounter);
checkInCodelets14054Ptr++;
(*checkInCodelets14055Ptr) = _checkInCodelets14055(1,1,this,codeletCounter);
checkInCodelets14055Ptr++;
(*checkInCodelets14052Ptr) = _checkInCodelets14052(1,1,this,codeletCounter);
checkInCodelets14052Ptr++;
(*checkInCodelets13997Ptr) = _checkInCodelets13997(1,1,this,codeletCounter);
(*checkInCodelets13997Ptr).decDep();
checkInCodelets13997Ptr++;
}
}
TP13995::~TP13995(){
delete [] TP14054_LoopCounterPerThread;delete [] TP14069_LoopCounterPerThread;delete [] barrierCodelets13995;
delete [] barrierCodelets14080;
delete [] checkInCodelets14080;
delete [] checkInCodelets14078;
delete [] barrierCodelets14077;
delete [] checkInCodelets14073;
delete [] checkInCodelets14069;
delete [] checkInCodelets14070;
delete [] checkInCodelets14067;
delete [] barrierCodelets14066;
delete [] checkInCodelets14058;
delete [] checkInCodelets14054;
delete [] checkInCodelets14055;
delete [] checkInCodelets14052;
delete [] barrierCodelets13997;
delete [] checkInCodelets13997;
}
/*TP13997: OMPForDirective*/
void TP13997::_barrierCodelets13997::fire(void)
{
TP13997* myTP = static_cast<TP13997*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets13997[0].decDep ();
}
bool TP13997::requestNewRangeIterations13997(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet13997 * codeletID;
int tempEndRange   = rangePerCodelet13997 * (codeletID + 1);
if (remainderRange13997 != 0)
{
if (codeletID < (uint32_t)remainderRange13997)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange13997;
tempEndRange += remainderRange13997;
}
}
tempStartRange = tempStartRange*1 + minIteration13997;
tempEndRange = tempEndRange*1 + minIteration13997;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration13997 < lastIteration13997)
{
(this->inputsTPParent->i_darts13997[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts13997[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration13997;
}
}
return isThereNewIteration;
}
void TP13997::_checkInCodelets13998::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 13998: ForStmt*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: m*/
int* i = &(this->inputsTPParent->i_darts13997[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts13997[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts13997[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts13997[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13997((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13997[0].decDep();
return;
}
for (int i_darts_counter_temp13997 = (*i);i_darts_counter_temp13997<=endRange && i_darts_counter_temp13997<=this->inputsTPParent->lastIteration13997;i_darts_counter_temp13997++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp13997 = (*j);
for(;j_darts_counter_temp13997 <= jend;j_darts_counter_temp13997++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp13997 = (*k);
for(;k_darts_counter_temp13997 <= nz - 2;k_darts_counter_temp13997++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp13997 = (*m);
for(;m_darts_counter_temp13997 < 5;m_darts_counter_temp13997++){
rsd[(i_darts_counter_temp13997)][j_darts_counter_temp13997][k_darts_counter_temp13997][m_darts_counter_temp13997] = dt * rsd[(i_darts_counter_temp13997)][j_darts_counter_temp13997][k_darts_counter_temp13997][m_darts_counter_temp13997];
}
(*m) = m_darts_counter_temp13997;
}
}
(*k) = k_darts_counter_temp13997;
}
}
(*j) = j_darts_counter_temp13997;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets13997[0].decDep();
}
TP13997::TP13997(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, int in_initIteration, int in_lastIteration, TP13997** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts13997(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts13997(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts13997(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts13997(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/, initIteration13997(in_initIteration), lastIteration13997(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets13997(new _barrierCodelets13997[1]) ,checkInCodelets13998(new _checkInCodelets13998[this->numThreads]){
/*Initialize the loop parameters*/
range13997 = abs (lastIteration13997 - initIteration13997) / 1;
rangePerCodelet13997 = range13997 / numThreads;
minIteration13997 = min<int>(lastIteration13997, initIteration13997);
remainderRange13997 = range13997 % numThreads;
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets13997[0] = _barrierCodelets13997(this->numThreads,this->numThreads,this, 0);
_checkInCodelets13998 * checkInCodelets13998Ptr = (this->checkInCodelets13998);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets13998);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets13998Ptr) = _checkInCodelets13998(2,1,this,codeletCounter);
#else
(*checkInCodelets13998Ptr) = _checkInCodelets13998(1,1,this,codeletCounter);
#endif
(*checkInCodelets13998Ptr).decDep();
checkInCodelets13998Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP13997::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets13998[localID].setID (codeletID);
this->checkInCodelets13998[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets13998[localID + this->baseNumThreads * i] = _checkInCodelets13998(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets13998[localID + this->baseNumThreads * i] = _checkInCodelets13998(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets13998[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets13998[localID + this->baseNumThreads * i].decDep();
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
TP13997::~TP13997(){
delete [] barrierCodelets13997;
delete [] checkInCodelets13998;
}
/*TP14054: ForStmt*/
void TP14054::_checkInCodelets14060::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif

/*printing node 14060: CallExpr*/
if (! __sync_val_compare_and_swap(&(myTP->controlTPParent->TP14060_alreadyLaunched), 0, 1))
{
/*Make the function call*/
invoke < TP_jacld > (myTP, myTP->numThreads, this->getID(), & (myTP->controlTPParent->checkInCodelets14061[this->getID()]) , & (myTP->controlTPParent->TPParent->barrierCodelets14066[0])  , &(myTP->controlTPParent->TP14060Ptr), (this->inputsTPParent->k_darts13995[this->getID()]));
}
else
{
if(myTP->controlTPParent->TP14060Ptr == nullptr)
{
myTP->add(this); 
return;
}
else
{
myTP->controlTPParent->TP14060Ptr->setNewInputs(
(this->inputsTPParent->k_darts13995[this->getID()]), this->getID());
myTP->controlTPParent->TP14060Ptr->nextCodeletsjacld[this->getID()] = &(myTP->controlTPParent->checkInCodelets14061[this->getID()]);
myTP->controlTPParent->TP14060Ptr->nextSyncCodeletsjacld[this->getID()] = &(myTP->controlTPParent->TPParent->barrierCodelets14066[0]);
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14060Ptr->firstCodelet[this->getID()].decDep();
#else
myTP->controlTPParent->TP14060Ptr->availableCodelets[this->getID()] = 1;
#endif
}
}
}
void TP14054::_checkInCodelets14061::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->k_darts14061= &(this->inputsTPParent->k_darts13995[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;

/*printing node 14062: CallExpr*/
printf("jacld complete... #k = %d\n", (*(this->inputsTPParent->k_darts14061)));
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling next codelet region: 14061 nextRegion: 14063 */
myTP->controlTPParent->checkInCodelets14063[this->getID()].decDep();
}
else
{
/*Signaling next codelet region: 14061 nextRegion: 14063 */
myTP->checkInCodelets14063[this->getID()].decDep();
}
}
void TP14054::_checkInCodelets14063::fire(void)
{

/*printing node 14063: CallExpr*/
if (! __sync_val_compare_and_swap(&(myTP->controlTPParent->TP14063_alreadyLaunched), 0, 1))
{
/*Make the function call*/
invoke < TP_blts > (myTP, myTP->numThreads, this->getID(), & (myTP->controlTPParent->checkInCodelets14064[this->getID()]) , & (myTP->controlTPParent->TPParent->barrierCodelets14066[0])  , &(myTP->controlTPParent->TP14063Ptr), nx, ny, nz, (this->inputsTPParent->k_darts13995[this->getID()]), omega, ist, iend, jst, jend, nx0, ny0);
}
else
{
if(myTP->controlTPParent->TP14063Ptr == nullptr)
{
myTP->add(this); 
return;
}
else
{
myTP->controlTPParent->TP14063Ptr->setNewInputs(
nx, ny, nz, (this->inputsTPParent->k_darts13995[this->getID()]), omega, ist, iend, jst, jend, nx0, ny0, this->getID());
myTP->controlTPParent->TP14063Ptr->nextCodeletsblts[this->getID()] = &(myTP->controlTPParent->checkInCodelets14064[this->getID()]);
myTP->controlTPParent->TP14063Ptr->nextSyncCodeletsblts[this->getID()] = &(myTP->controlTPParent->TPParent->barrierCodelets14066[0]);
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14063Ptr->firstCodelet[this->getID()].decDep();
#else
myTP->controlTPParent->TP14063Ptr->availableCodelets[this->getID()] = 1;
#endif
}
}
}
void TP14054::_checkInCodelets14064::fire(void)
{
if (this->getID() == 0)
{
/*Init the vars for this region*/
/*Initialize the vars of the inlined region*/
this->inputsTPParent->k_darts14064= &(this->inputsTPParent->k_darts13995[this->getLocalID()])/*OMP_SHARED_PRIVATE - VAR INLINED*/;

/*printing node 14065: CallExpr*/
printf("blts complete... #k = %d\n", (*(this->inputsTPParent->k_darts14064)));
/*Signaling next codelet from last stmt in the codelet*/
/*The node is the last one in a complex loop, so signal the inc node*/
myTP->controlTPParent->TPParent->checkInCodelets14058[this->getID()].decDep();
}
else
{
/*The node is the last one in a complex loop, so signal the inc node*/
myTP->TPParent->checkInCodelets14058[this->getID()].decDep();
}
}
TP14054::TP14054(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, TP13995* in_inputsTPParent, TP14054** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(in_inputsTPParent), ptrToThisTP(in_ptrToThisTP), TP14060Ptr(nullptr), TP14060_alreadyLaunched(0), TP14061_alreadyLaunched(0), TP14063Ptr(nullptr), TP14063_alreadyLaunched(0), TP14064_alreadyLaunched(0) ,checkInCodelets14060(new _checkInCodelets14060[this->numThreads]) ,checkInCodelets14061(new _checkInCodelets14061[this->numThreads]) ,checkInCodelets14063(new _checkInCodelets14063[this->numThreads]) ,checkInCodelets14064(new _checkInCodelets14064[this->numThreads]){
/*Initialize Codelets*/
_checkInCodelets14064 * checkInCodelets14064Ptr = (this->checkInCodelets14064);
_checkInCodelets14063 * checkInCodelets14063Ptr = (this->checkInCodelets14063);
_checkInCodelets14061 * checkInCodelets14061Ptr = (this->checkInCodelets14061);
_checkInCodelets14060 * checkInCodelets14060Ptr = (this->checkInCodelets14060);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets14060);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets14064Ptr) = _checkInCodelets14064(1,1,this,codeletCounter);
checkInCodelets14064Ptr++;
(*checkInCodelets14063Ptr) = _checkInCodelets14063(1,1,this,codeletCounter);
checkInCodelets14063Ptr++;
(*checkInCodelets14061Ptr) = _checkInCodelets14061(1,1,this,codeletCounter);
checkInCodelets14061Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets14060Ptr) = _checkInCodelets14060(2,1,this,codeletCounter);
#else
(*checkInCodelets14060Ptr) = _checkInCodelets14060(1,1,this,codeletCounter);
#endif
(*checkInCodelets14060Ptr).decDep();
checkInCodelets14060Ptr++;
}
*(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
this->firstCodelet[this->getID()].decDep();
#else
this->availableCodelets[this->getID()] = 1;
#endif
}
TP14054::~TP14054(){
delete [] checkInCodelets14064;
delete [] checkInCodelets14063;
delete [] checkInCodelets14061;
delete [] checkInCodelets14060;
}
/*TP14069: ForStmt*/
void TP14069::_checkInCodelets14075::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif

/*printing node 14075: CallExpr*/
if (! __sync_val_compare_and_swap(&(myTP->controlTPParent->TP14075_alreadyLaunched), 0, 1))
{
/*Make the function call*/
invoke < TP_jacu > (myTP, myTP->numThreads, this->getID(), & (myTP->controlTPParent->checkInCodelets14076[this->getID()]) , & (myTP->controlTPParent->TPParent->barrierCodelets14077[0])  , &(myTP->controlTPParent->TP14075Ptr), (this->inputsTPParent->k_darts13995[this->getID()]));
}
else
{
if(myTP->controlTPParent->TP14075Ptr == nullptr)
{
myTP->add(this); 
return;
}
else
{
myTP->controlTPParent->TP14075Ptr->setNewInputs(
(this->inputsTPParent->k_darts13995[this->getID()]), this->getID());
myTP->controlTPParent->TP14075Ptr->nextCodeletsjacu[this->getID()] = &(myTP->controlTPParent->checkInCodelets14076[this->getID()]);
myTP->controlTPParent->TP14075Ptr->nextSyncCodeletsjacu[this->getID()] = &(myTP->controlTPParent->TPParent->barrierCodelets14077[0]);
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14075Ptr->firstCodelet[this->getID()].decDep();
#else
myTP->controlTPParent->TP14075Ptr->availableCodelets[this->getID()] = 1;
#endif
}
}
}
void TP14069::_checkInCodelets14076::fire(void)
{

/*printing node 14076: CallExpr*/
if (! __sync_val_compare_and_swap(&(myTP->controlTPParent->TP14076_alreadyLaunched), 0, 1))
{
/*Make the function call*/
invoke < TP_buts > (myTP, myTP->numThreads, this->getID(), & (myTP->controlTPParent->TPParent->checkInCodelets14073[this->getID()]) , & (myTP->controlTPParent->TPParent->barrierCodelets14077[0])  , &(myTP->controlTPParent->TP14076Ptr), nx, ny, nz, (this->inputsTPParent->k_darts13995[this->getID()]), omega, ist, iend, jst, jend, nx0, ny0);
}
else
{
if(myTP->controlTPParent->TP14076Ptr == nullptr)
{
myTP->add(this); 
return;
}
else
{
myTP->controlTPParent->TP14076Ptr->setNewInputs(
nx, ny, nz, (this->inputsTPParent->k_darts13995[this->getID()]), omega, ist, iend, jst, jend, nx0, ny0, this->getID());
myTP->controlTPParent->TP14076Ptr->nextCodeletsbuts[this->getID()] = &(myTP->controlTPParent->TPParent->checkInCodelets14073[this->getID()]);
myTP->controlTPParent->TP14076Ptr->nextSyncCodeletsbuts[this->getID()] = &(myTP->controlTPParent->TPParent->barrierCodelets14077[0]);
#if USE_SPIN_CODELETS == 0
myTP->controlTPParent->TP14076Ptr->firstCodelet[this->getID()].decDep();
#else
myTP->controlTPParent->TP14076Ptr->availableCodelets[this->getID()] = 1;
#endif
}
}
}
TP14069::TP14069(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, TP13995* in_inputsTPParent, TP14069** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(in_inputsTPParent), ptrToThisTP(in_ptrToThisTP), TP14075Ptr(nullptr), TP14075_alreadyLaunched(0), TP14076Ptr(nullptr), TP14076_alreadyLaunched(0) ,checkInCodelets14075(new _checkInCodelets14075[this->numThreads]) ,checkInCodelets14076(new _checkInCodelets14076[this->numThreads]){
/*Initialize Codelets*/
_checkInCodelets14076 * checkInCodelets14076Ptr = (this->checkInCodelets14076);
_checkInCodelets14075 * checkInCodelets14075Ptr = (this->checkInCodelets14075);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets14075);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++)
{
(*checkInCodelets14076Ptr) = _checkInCodelets14076(1,1,this,codeletCounter);
checkInCodelets14076Ptr++;
#if USE_SPIN_CODELETS == 0
(*checkInCodelets14075Ptr) = _checkInCodelets14075(2,1,this,codeletCounter);
#else
(*checkInCodelets14075Ptr) = _checkInCodelets14075(1,1,this,codeletCounter);
#endif
(*checkInCodelets14075Ptr).decDep();
checkInCodelets14075Ptr++;
}
*(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
this->firstCodelet[this->getID()].decDep();
#else
this->availableCodelets[this->getID()] = 1;
#endif
}
TP14069::~TP14069(){
delete [] checkInCodelets14076;
delete [] checkInCodelets14075;
}
/*TP14080: OMPForDirective*/
void TP14080::_barrierCodelets14080::fire(void)
{
TP14080* myTP = static_cast<TP14080*>(myTP_);
myTP->controlTPParent->TPParent->barrierCodelets14080[0].decDep ();
}
bool TP14080::requestNewRangeIterations14080(int* endRange, uint32_t codeletID)
{
/*Scheduling Policy = Static */
/*Chunk = 0*/
bool isThereNewIteration = false;
{
/*Static Scheduling*/
int tempStartRange = rangePerCodelet14080 * codeletID;
int tempEndRange   = rangePerCodelet14080 * (codeletID + 1);
if (remainderRange14080 != 0)
{
if (codeletID < (uint32_t)remainderRange14080)
{
tempStartRange += codeletID;
tempEndRange += (codeletID + 1);
}
else
{
tempStartRange += remainderRange14080;
tempEndRange += remainderRange14080;
}
}
tempStartRange = tempStartRange*1 + minIteration14080;
tempEndRange = tempEndRange*1 + minIteration14080;
if(tempStartRange != tempEndRange)
{
isThereNewIteration = true;
}
if(initIteration14080 < lastIteration14080)
{
(this->inputsTPParent->i_darts14080[codeletID]) = min(tempStartRange, tempEndRange);
*endRange   = max(tempStartRange, tempEndRange);
}
else
{
(this->inputsTPParent->i_darts14080[codeletID]) = max(tempStartRange, tempEndRange);
*endRange   = min(tempStartRange, tempEndRange);
}
if(codeletID == this->numThreads - 1)
{
*endRange = *endRange + 1;
}
if(codeletID == this->numThreads - 1)
{
*endRange = lastIteration14080;
}
}
return isThereNewIteration;
}
void TP14080::_checkInCodelets14081::fire(void)
{
#if USE_SPIN_CODELETS == 1
/*Wait until the codelet with the same ID finishes in the previous TP*/
if(myTP->availableCodelets[this->getLocalID()] == 0)
{
myTP->add(this);
return;
}
#endif
/*Init the vars for this region*/

/*printing node 14081: ForStmt*/
/*var: i*/
/*var: j*/
/*var: k*/
/*var: m*/
/*var: tmp*/
int* i = &(this->inputsTPParent->i_darts14080[this->getLocalID()]);
(void)i/*OMP_PRIVATE*/;
int* j = &(this->inputsTPParent->j_darts14080[this->getLocalID()]);
(void)j/*OMP_PRIVATE*/;
int* k = &(this->inputsTPParent->k_darts14080[this->getLocalID()]);
(void)k/*OMP_PRIVATE*/;
int* m = &(this->inputsTPParent->m_darts14080[this->getLocalID()]);
(void)m/*OMP_PRIVATE*/;
double* tmp = (this->inputsTPParent->tmp_darts14080);
(void)tmp/*OMP_SHARED*/;
bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14080((int*)&(this->endRange), this->getLocalID());if(isThereNewIteration == false)
{
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets14080[0].decDep();
return;
}
for (int i_darts_counter_temp14080 = (*i);i_darts_counter_temp14080<=endRange && i_darts_counter_temp14080<=this->inputsTPParent->lastIteration14080;i_darts_counter_temp14080++)
{
{
{
/*Loop's init*/
(*j) = jst;
int j_darts_counter_temp14080 = (*j);
for(;j_darts_counter_temp14080 <= jend;j_darts_counter_temp14080++){
{
/*Loop's init*/
(*k) = 1;
int k_darts_counter_temp14080 = (*k);
for(;k_darts_counter_temp14080 <= nz - 2;k_darts_counter_temp14080++){
{
/*Loop's init*/
(*m) = 0;
int m_darts_counter_temp14080 = (*m);
for(;m_darts_counter_temp14080 < 5;m_darts_counter_temp14080++){
u[(i_darts_counter_temp14080)][j_darts_counter_temp14080][k_darts_counter_temp14080][m_darts_counter_temp14080] = u[(i_darts_counter_temp14080)][j_darts_counter_temp14080][k_darts_counter_temp14080][m_darts_counter_temp14080] + (*(tmp)) * rsd[(i_darts_counter_temp14080)][j_darts_counter_temp14080][k_darts_counter_temp14080][m_darts_counter_temp14080];
}
(*m) = m_darts_counter_temp14080;
}
}
(*k) = k_darts_counter_temp14080;
}
}
(*j) = j_darts_counter_temp14080;
}
}
}
/*Signaling next codelet from last stmt in the codelet*/
/*Signaling omp for stmt's barrier*/
myTP->controlTPParent->barrierCodelets14080[0].decDep();
}
TP14080::TP14080(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, int in_initIteration, int in_lastIteration, double *in_tmp, TP14080** in_ptrToThisTP):ompTP(in_numThreads, in_mainCodeletID), TPParent(in_TPParent), controlTPParent(this), inputsTPParent(this),i_darts14080(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,j_darts14080(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,k_darts14080(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,m_darts14080(new int[this->numThreads])/*OMP_PRIVATE - INPUT*/,tmp_darts14080(in_tmp)/*OMP_SHARED - INPUT*/, initIteration14080(in_initIteration), lastIteration14080(in_lastIteration), readyCodelets(0), baseNumThreads(this->numThreads / DARTS_CODELETS_MULT) ,barrierCodelets14080(new _barrierCodelets14080[1]) ,checkInCodelets14081(new _checkInCodelets14081[this->numThreads]){
/*Initialize the loop parameters*/
range14080 = abs (lastIteration14080 - initIteration14080) / 1;
rangePerCodelet14080 = range14080 / numThreads;
minIteration14080 = min<int>(lastIteration14080, initIteration14080);
remainderRange14080 = range14080 % numThreads;
/*Initialize inputs and vars.*/
/*Initialize Codelets*/
barrierCodelets14080[0] = _barrierCodelets14080(this->numThreads,this->numThreads,this, 0);
_checkInCodelets14081 * checkInCodelets14081Ptr = (this->checkInCodelets14081);
#if USE_SPIN_CODELETS == 0
firstCodelet = (this->checkInCodelets14081);
#endif
for(size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads; codeletCounter++)
{
#if USE_SPIN_CODELETS == 0
(*checkInCodelets14081Ptr) = _checkInCodelets14081(2,1,this,codeletCounter);
#else
(*checkInCodelets14081Ptr) = _checkInCodelets14081(1,1,this,codeletCounter);
#endif
(*checkInCodelets14081Ptr).decDep();
checkInCodelets14081Ptr++;
}
this->dispatchCodelet(this->getID());
*(in_ptrToThisTP) = this;
}
void TP14080::dispatchCodelet(size_t codeletID)
{
int idx = codeletID / this->baseNumThreads;
int localID = codeletID - this->baseNumThreads * idx;
this->checkInCodelets14081[localID].setID (codeletID);
this->checkInCodelets14081[localID].setLocalID (localID);
/*Check if we want to replicate codelets*/
for(size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++){
#if USE_SPIN_CODELETS == 0
this->checkInCodelets14081[localID + this->baseNumThreads * i] = _checkInCodelets14081(2, 1, this, localID + this->baseNumThreads * i);
#else
this->checkInCodelets14081[localID + this->baseNumThreads * i] = _checkInCodelets14081(1, 1, this, localID + this->baseNumThreads * i);
#endif
 this->checkInCodelets14081[localID + this->baseNumThreads * i].setID(codeletID);
 this->checkInCodelets14081[localID + this->baseNumThreads * i].decDep();
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
TP14080::~TP14080(){
delete [] barrierCodelets14080;
delete [] checkInCodelets14081;
}
