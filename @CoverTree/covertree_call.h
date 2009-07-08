#ifndef COVERTREE_CALL_H
#define COVERTREE_CALL_H

#include <mex.h>

struct CallbackParam;

static double distanceMlevalCallback(mxArray *p1, mxArray *p2,
                                     double upper_bound);
static double distanceMlfuncCallback(mxArray *p1, mxArray *p2,
                                     double upper_bound);
static void serializeMlevalCallback(int nlhs, mxArray *plhs[],
                                    int nrhs, const mxArray *prhs[],
                                    const CallbackParam &param);
static void serializeMlfuncCallback(int nlhs, mxArray *plhs[],
                                    int nrhs, const mxArray *prhs[],
                                    const CallbackParam &param);

static void ctBatchCreate(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[]);
static void ctInsert(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[]);
static void ctKNearestNeighbor(int nlhs, mxArray *plhs[],
                               int nrhs, const mxArray *prhs[]);
static void ctEpsilonNearestNeighbor(int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[]);
static void ctUnequalNearestNeighbor(const int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[]);
static void ctDepthDist(int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[]);
static void ctHeightDist(int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[]);
static void ctBreadthDist(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[]);
static void ctDelete(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[]);
static void ctLoad(int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[]);
static void ctSave(int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[]);

#endif
