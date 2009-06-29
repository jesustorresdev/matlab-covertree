// distances.cc - Some distance functions to use with CoverTree

#include <mex.h>
#include <stdint.h>
#include <values.h>

#include "distances.h"

const int BATCH = 152;

// distance:  euclidean
// type:      real vectors
// precision: double
// labelled:  false
extern "C" float dist_ervdu(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<double, 0, BATCH>(p1, p2, upper_bound);
}

// distance:  euclidean
// type:      real vectors
// precision: double
// labelled:  true
extern "C" float dist_ervdl(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<double, 1, BATCH>(p1, p2, upper_bound);
}

// distance:  euclidean
// type:      real vectors
// precision: single
// labelled:  false
extern "C" float dist_ervsu(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<float, 0, BATCH>(p1, p2, upper_bound);
}

// distance:  euclidean
// type:      real vectors
// precision: single
// labelled:  true
extern "C" float dist_ervsl(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<float, 1, BATCH>(p1, p2, upper_bound);
}

// distance:  euclidean
// type:      real vectors
// precision: uint8
// labelled:  false
extern "C" float dist_ervbu(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<uint8_t, 0, BATCH>(p1, p2, upper_bound);
}

// distance:  euclidean
// type:      real vectors
// precision: uint8
// labelled:  true
extern "C" float dist_ervbl(const mxArray *p1,
                            const mxArray *p2,
                            float upper_bound)
{
  return vectorsEuclideanDistance<uint8_t, 1, BATCH>(p1, p2, upper_bound);
}

// Entry point to C MEX
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mexWarnMsgTxt("This MEX is only to use through CoverTree objects.");

  if (nrhs < 2)
    mexErrMsgTxt("Too few input arguments.");

  const mxArray *p1 = prhs[0], *p2 = prhs[1];
  float d, upper_bound = MAXFLOAT;

  if (nrhs > 2)
      upper_bound = mxGetScalar(prhs[2]);

  switch(mxGetClassID(p1)) {
    case mxDOUBLE_CLASS:
      d = dist_ervdu(p1, p2, upper_bound);
      break;

    case mxSINGLE_CLASS:
      d = dist_ervsu(p1, p2, upper_bound);
      break;

    case mxUINT8_CLASS:
      d = dist_ervbu(p1, p2, upper_bound);
      break;

    default:
      plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
      return;
  }

  plhs[0] = mxCreateDoubleScalar(d);
}
