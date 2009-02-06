// distances.cc - Some distance functions to use with CoverTree

#include <alloca.h>
#include <blas.h>
#include <math.h>
#include <mex.h>
#include <stdint.h>
#include <values.h>

// Basic vector-vector operations
//

// MINUS: z <- x-y
template<class T>
inline
__attribute__((optimize("tree-vectorize")))
void minus(T * __restrict__ vd, const T *v1, const T *v2, int n)
{
  for (int i = 0; i < n; i++)
    vd[i] = v1[i] - v2[i];
}

// DOT: dot <- x'y
template<class T>
__attribute__((optimize("tree-vectorize")))
inline
double dot(const T *v1, const T *v2, int n)
{
  double sum = 0.;
  for (int i = 0; i < n; i++)
    sum += v1[i] * v2[i];
  return sum;
}

template<>
inline
double dot(const float *v1, const float *v2, int n)
{
  int inc = 1;
  return sdot_(&n, (float*)v1, &inc, (float*)v2, &inc);
}

template<>
inline
double dot(const double *v1, const double *v2, int n)
{
  int inc = 1;
  return ddot_(&n, (double*)v1, &inc, (double*)v2, &inc);
}

// Euclidean distances
//

template<class T, int offset, int batch>
inline
float vectorsEuclideanDistance(const mxArray *pm_p1,
                               const mxArray *pm_p2,
                               float upper_bound)
{
  const T *p2 = (const T*)mxGetData(pm_p2) + offset;
  const T *p1 = (const T*)mxGetData(pm_p1) + offset;
  int len = mxGetNumberOfElements(pm_p1) - offset;

  double sum = 0;
  T * __restrict__ tmp = (T*)alloca(batch * sizeof(T));

  const T *end = p1 + len;
  upper_bound *= upper_bound;
  for (const T *batch_end = p1 + batch; batch_end <= end; batch_end += batch) {
    minus<T>(tmp, p1, p2, batch);
    sum += dot<T>(tmp, tmp, batch);
    if (sum > upper_bound)
      return sqrt(sum);

    p1 = batch_end;
    p2 += batch; len -= batch;
  }

  minus<T>(tmp, p1, p2, len);
  sum += dot<T>(tmp, tmp, len);
  return sqrt(sum);
}

// Exported distance functions
//

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
