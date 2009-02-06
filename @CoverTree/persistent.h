#ifndef _PERSISTENT_H
#define _PERSISTENT_H
#include <mex.h>
#include "stack.h"

namespace persistent {

v_array<mxArray*> arrays;
v_array<node<mxArray*>* > children;

template<class T>
void make_persistent(v_array<T> &v)
{
  T *ele_end = v.elements + v.index;
  for (T *ele = v.elements; ele != ele_end; ele++)
    internal_make_persistent(*ele);
  v.index = 0;
}

template<class T>
inline
void internal_make_persistent(T *pr)
{
  mexMakeMemoryPersistent(pr);
}

template<>
inline
void internal_make_persistent(mxArray *pm)
{
  mexMakeArrayPersistent(pm);
}

template<class T>
void free(v_array<T> &v)
{
  T *ele_end = v.elements + v.index;
  for (T *ele = v.elements; ele != ele_end; ele++)
    internal_free(*ele);
  ::free(v);
}

template<class T>
inline
void internal_free(T *pr)
{
  mxFree(pr);
}

template<>
inline
void internal_free(mxArray *pm)
{
  mxDestroyArray(pm);
}

}

#endif
