// mex_undoc.h - Undocumented MEX functions

#ifndef MEX_SERIALIZATION_H
#define MEX_SERIALIZATION_H

#include <mex.h>

EXTERN_C mxArray* mxCreateReference(const mxArray *pm);

// Serialization / Deserialization
EXTERN_C mxArray* mxSerialize(const mxArray *pm);
EXTERN_C mxArray* mxDeserialize(const void *pr, mwSize n);

#endif
