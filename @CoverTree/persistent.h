// persistent.h
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
