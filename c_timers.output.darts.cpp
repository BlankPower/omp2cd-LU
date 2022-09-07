#include "c_timers.output.darts.h"
using namespace darts;
using namespace std;
double elapsed_dartsThreadPriv0[MAXNUMTHREADS][64];
double start_dartsThreadPriv0[MAXNUMTHREADS][64];
/*Function: elapsed_time, ID: 39*/
double elapsed_time()
{
    /*elapsed_time:39*/
    /*CompoundStmt:176*/
    double t;
    wtime(&t);
    return (t);
}
/*Function: timer_clear, ID: 40*/
void timer_clear(int n)
{
    /*timer_clear:40*/
    /*CompoundStmt:181*/
    ((elapsed_dartsThreadPriv0[0]))[n] = 0.;
}
/*Function: timer_start, ID: 41*/
void timer_start(int n)
{
    /*timer_start:41*/
    /*CompoundStmt:185*/
    ((start_dartsThreadPriv0[0]))[n] = elapsed_time();
}
/*Function: timer_stop, ID: 42*/
void timer_stop(int n)
{
    /*timer_stop:42*/
    /*CompoundStmt:189*/
    double t, now;
    now = elapsed_time();
    t = now - ((start_dartsThreadPriv0[0]))[n];
    ((elapsed_dartsThreadPriv0[0]))[n] += t;
}
/*Function: timer_read, ID: 43*/
double timer_read(int n)
{
    /*timer_read:43*/
    /*CompoundStmt:198*/
    return (((elapsed_dartsThreadPriv0[0]))[n]);
}
