#include "c_print_results.output.darts.h"
using namespace darts;
using namespace std;
/*Function: c_print_results, ID: 37*/
void c_print_results(char* name, char class_is, int n1, int n2, int n3, int niter, int nthreads,
    double t, double mops, char* optype, int passed_verification, char* npbversion,
    char* compiletime, char* cc, char* clink, char* c_lib, char* c_inc, char* cflags,
    char* clinkflags, char* rand)
{
    /*c_print_results:37*/
    /*CompoundStmt:114*/
    char* evalue = "1000";
    printf("\n\n %s Benchmark Completed\n", name);
    printf(" class_is           =                        %c\n", class_is);
    if (n2 == 0 && n3 == 0)
        printf(" Size            =             %12d\n", n1);
    else
        printf(" Size            =              %3dx%3dx%3d\n", n1, n2, n3);
    printf(" Iterations      =             %12d\n", niter);
    printf(" Threads         =             %12d\n", nthreads);
    printf(" Time in seconds =             %12.2f\n", t);
    printf(" Mop/s total     =             %12.2f\n", mops);
    printf(" Operation type  = %24s\n", optype);
    if (passed_verification)
        printf(" Verification    =               SUCCESSFUL\n");
    else
        printf(" Verification    =             UNSUCCESSFUL\n");
    printf(" Version         =           %12s\n", npbversion);
    printf(" Compile date    =             %12s\n", compiletime);
    printf("\n Compile options:\n");
    printf("    CC           = %s\n", cc);
    printf("    CLINK        = %s\n", clink);
    printf("    C_LIB        = %s\n", c_lib);
    printf("    C_INC        = %s\n", c_inc);
    printf("    CFLAGS       = %s\n", cflags);
    printf("    CLINKFLAGS   = %s\n", clinkflags);
    printf("    RAND         = %s\n", rand);
}
