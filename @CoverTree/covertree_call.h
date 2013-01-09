// covertree_call.h
//
//   Written by Jesús Torres <jmtorres@ull.es>
//
// This file is part of MATLAB Binding for Cover Tree.
//
// MATLAB Binding for Cover Tree is free software: you can redistribute it
// and/or modify it under the terms of the GNU General Public License or
// Lesser GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// MATLAB Binding for Cover Tree is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// (Lesser) GNU General Public License for more details.
//
// You should have received a copy of the (Lesser) GNU General Public License
// along with MATLAB Binding for Cover Tree.  If not, see
// <http://www.gnu.org/licenses/>.
//

#ifndef COVERTREE_CALL_H
#define COVERTREE_CALL_H

#include <mex.h>

struct CallbackParam;

static float distanceMlevalCallback(mxArray *p1, mxArray *p2,
                                    float upper_bound);
static float distanceMlfuncCallback(mxArray *p1, mxArray *p2,
                                    float upper_bound);
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
