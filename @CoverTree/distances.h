// distances.h - Some helper functions to distance computation

#include <alloca.h>
#include <blas.h>
#include <math.h>

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
