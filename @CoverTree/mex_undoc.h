// mex_undoc.h - Undocumented MEX functions
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

#ifndef MEX_SERIALIZATION_H
#define MEX_SERIALIZATION_H

#include <mex.h>

EXTERN_C mxArray* mxCreateReference(const mxArray *pm);

// Serialization / Deserialization
EXTERN_C mxArray* mxSerialize(const mxArray *pm);
EXTERN_C mxArray* mxDeserialize(const void *pr, mwSize n);

#endif
