// stack.h
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

#ifndef _STACK_H
#define _STACK_H
#include <mex.h>
#include <stdlib.h>

template<class T> class v_array{
 public:
  int index;
  int length;
  T* elements;

  T last() { return elements[index-1];}
  void decr() { index--;}
  v_array() { index = 0; length = 0; elements = NULL;}
  T& operator[](unsigned int i) { return elements[i]; }
};

template<class T> void push(v_array<T>& v, const T &new_ele)
{
  while(v.index >= v.length)
    {
      v.length = 2*v.length + 3;
      v.elements = (T *)mxRealloc(v.elements,sizeof(T) * v.length);
    }
  v[v.index++] = new_ele;
}

template<class T> void alloc(v_array<T>& v, int length)
{
  v.elements = (T *)mxRealloc(v.elements, sizeof(T) * length);
  v.length = length;
}

template<class T> v_array<T> pop(v_array<v_array<T> > &stack)
{
  if (stack.index > 0)
    return stack[--stack.index];
  else
    return v_array<T>();
}

template<class T> void free(v_array<T>& v)
{
  mxFree(v.elements);
  v.index = 0; v.length=0; v.elements = NULL;
}

#endif
