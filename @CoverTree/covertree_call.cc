// covertree_call.cc - MEX wrapper of Cover Tree functions

#include <mex.h>

#include <dlfcn.h>
#include <stdint.h>
#include <string.h>

#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/split_free.hpp>

#include "cover_tree.h"
#include "covertree_call.h"
#include "mex_undoc.h"
#include "persistent.h"
#include "stack.h"

#define COVER_TREE_CLASS               "CoverTree"
const char PROP_COVER_TREE_HANDLER[] = "CoverTreeHandle";
const char OBJECTNAME_DELIMITER      = ':';

// CallbackParam store the parameters for a callback function
struct CallbackParam
{
  mxArray *mleval_function;
  const char *mlfunc_function;
  void *object_handle;

  CallbackParam() :
      mleval_function(NULL), mlfunc_function(NULL), object_handle(NULL) {}
};

typedef static double (*DistanceCallback)(mxArray*, mxArray*, double);
typedef static void (*MexCallback)(int, mxArray*[], int, const mxArray*[],
                                   const CallbackParam&);

// Properties of CoverTree MATLAB object that containts a callback function
enum CoverTreeFcnProperties {
  PROP_DISTANCE_FCN = 0,
  PROP_PRESERIALIZE_FCN = 1,
  PROP_POSTDESERIALIZE_FCN = 2,
  PROP_END = 2,
};

// Pending dynamic library handles
void *PendingHandles[] = {
  NULL,                         // PROP_DISTANCE_FCN
  NULL,                         // PROP_PRESERIALIZE_FCN
  NULL,                         // PROP_POSTDESERIALIZE_FCN
};

struct CoverTree
{
  const mxArray* pm_ctobj;
  node<mxArray*> root_node;

  template <class C>
  struct Callback
  {
    C function;
    CallbackParam param;
    CoverTreeFcnProperties property;

    Callback(CoverTreeFcnProperties prop) : property(prop) {}
  };

  CoverTree::Callback<DistanceCallback> distance;
  CoverTree::Callback<MexCallback> pre_serialize;
  CoverTree::Callback<MexCallback> post_deserialize;
};

// Current CoverTree must be global to allow access to callback functions
static CoverTree *CT;
// Call the callback function DistanceFcn when distance() is called
static DistanceCallback distance = NULL;

struct PropertyCallbacks
{
  const char *name;

  union Callback{
    DistanceCallback distance;
    MexCallback mex;

    Callback(DistanceCallback cb) : distance(cb) {}
    Callback(MexCallback cb) : mex(cb) {}

    Callback& operator =(DistanceCallback cb) {distance = cb; return *this;}
    Callback& operator =(MexCallback cb) {mex = cb; return *this;}

    operator DistanceCallback() {return distance;}
    operator MexCallback() {return mex;}
  };

  PropertyCallbacks::Callback mleval_callback;
  PropertyCallbacks::Callback mlfunc_callback;
};

static PropertyCallbacks callbacks[] = {
  {"DistanceFcn", distanceMlevalCallback, distanceMlfuncCallback},
  {"PreSerializeFcn", serializeMlevalCallback, serializeMlfuncCallback},
  {"PostDeserializeFcn", serializeMlevalCallback, serializeMlfuncCallback},
};

// Cover Tree wrappered functions
struct Function
{
  const char *name;
  void (*func)(int, mxArray *[], int, const mxArray *[]);
};

static Function functions[] = {
  {"batch_create", ctBatchCreate},
  {"insert", ctInsert},
  {"k_nearest_neighbor", ctKNearestNeighbor},
  {"epsilon_nearest_neighbor", ctEpsilonNearestNeighbor},
  {"unequal_nearest_neighbor", ctUnequalNearestNeighbor},
  {"depth_dist", ctDepthDist},
  {"height_dist", ctHeightDist},
  {"breadth_dist", ctBreadthDist},
  {"delete", ctDelete},
  {"load", ctLoad},
  {"save", ctSave},
  {NULL, NULL}
};

static double distanceMlevalCallback(mxArray *p1, mxArray *p2,
                                     double upper_bound)
{
  mxArray *plhs = NULL;
  mxArray *ub = mxCreateDoubleScalar(upper_bound);
  mxArray *prhs[4] = {CT->distance.param.mleval_function, p1, p2, ub};
  mexCallMATLAB(1, &plhs, 4, prhs, "feval");

  if (! mxIsNumeric(plhs) || mxGetNumberOfElements(plhs) != 1) {
    mexErrMsgIdAndTxt(COVER_TREE_CLASS":distance:ReturnError",
                      "Distance function must return a numeric scalar");
  }

  double dist = mxGetScalar(plhs);
  mxDestroyArray(ub);
  mxDestroyArray(plhs);

  return dist;
}

static double distanceMlfuncCallback(mxArray *p1, mxArray *p2,
                                     double upper_bound)
{
  mxArray *plhs = NULL;
  mxArray *ub = mxCreateDoubleScalar(upper_bound);
  mxArray *prhs[3] = {p1, p2, ub};
  mexCallMATLAB(1, &plhs, 3, prhs, CT->distance.param.mlfunc_function);

  if (! mxIsNumeric(plhs) || mxGetNumberOfElements(plhs) != 1) {
    mexErrMsgIdAndTxt(COVER_TREE_CLASS":distance:ReturnError",
                      "Distance function '%s' must return a numeric scalar",
                      CT->distance.param.mlfunc_function);
  }

  double dist = mxGetScalar(plhs);
  mxDestroyArray(ub);
  mxDestroyArray(plhs);

  return dist;
}

static void serializeMlevalCallback(int nlhs, mxArray *plhs[],
                                    int nrhs, const mxArray *prhs[],
                                    const CallbackParam &param)
{
  mxArray *new_prhs[++nrhs];
  new_prhs[0] = param.mleval_function;
  memcpy(new_prhs + 1, prhs, sizeof(new_prhs));
  mexCallMATLAB(nlhs, plhs, nrhs, new_prhs, "feval");
}

static void serializeMlfuncCallback(int nlhs, mxArray *plhs[],
                                    int nrhs, const mxArray *prhs[],
                                    const CallbackParam &param)
{
  mexCallMATLAB(nlhs, plhs, nrhs, (mxArray**)prhs, param.mlfunc_function);
}

template<class C>
static CoverTree::Callback<C> getCoverTreeCallback(const mxArray *pm,
                                                   CoverTreeFcnProperties prop)
{
  CoverTree::Callback<C> callback(prop);
  mxArray *pv = mxGetProperty(pm, 0, callbacks[prop].name);

  if (mxIsEmpty(pv))
    callback.function = NULL;
  else
  {
    if (mxGetClassID(pv) == mxFUNCTION_CLASS) {
      callback.function = callbacks[prop].mleval_callback;
      callback.param.mleval_function = pv;
    }
    else if (mxGetClassID(pv) == mxCHAR_CLASS) {
      char *s1 =  mxArrayToString(pv);
      char *s2 = strchr(s1, OBJECTNAME_DELIMITER);
      if (s2 == NULL) {
        callback.function = callbacks[prop].mlfunc_callback;
        callback.param.mlfunc_function = s1;
        mxDestroyArray(pv);
      }
      else {
        *(s2++) = '\0';
        void *handle = dlopen(s1, RTLD_NOW | RTLD_LOCAL);
        if(handle == NULL)
          mexErrMsgIdAndTxt(COVER_TREE_CLASS":dlopen:error", dlerror());

        dlerror();
        void *sym = dlsym(handle, s2);
        char *error = dlerror();
        if(error != NULL)
          mexErrMsgIdAndTxt(COVER_TREE_CLASS":dlsym:error", error);

        callback.function = (C)sym;
        callback.param.object_handle = PendingHandles[prop] = handle;

        mxFree(s1);
        mxDestroyArray(pv);
      }
    }
    else {
      mexErrMsgIdAndTxt(COVER_TREE_CLASS":InvalidCallback",
                        "Invalid callback function in property '%s'.",
                        callbacks[prop].name);
    }
  }

  return callback;
}

template<class C>
static void freeCoverTreeCallback(CoverTree::Callback<C> &callback)
{
  if (callback.param.mleval_function)
    mxDestroyArray(callback.param.mleval_function);
  else if (callback.param.mlfunc_function)
    mxFree((void*)callback.param.mlfunc_function);
  else if(callback.param.object_handle != NULL)
    dlclose(callback.param.object_handle);
}

template<class C>
static void makeCoverTreeCallbackPersistent(CoverTree::Callback<C> &callback)
{
  if (callback.param.mleval_function)
    mexMakeArrayPersistent(callback.param.mleval_function);
  else if (callback.param.mlfunc_function)
    mexMakeMemoryPersistent((void*)callback.param.mlfunc_function);
  else if(callback.param.object_handle != NULL)
    PendingHandles[callback.property] = NULL;
}

static CoverTree* getCoverTree(const mxArray *pm)
{
  if (! mxIsClass(pm, COVER_TREE_CLASS))
    mexErrMsgTxt("A "COVER_TREE_CLASS" object is required.");

  CoverTree *ct = NULL;

  // Get the CoverTree struct pointer
  mxArray *pm_ct = mxGetProperty(pm, 0, PROP_COVER_TREE_HANDLER);
  if (mxIsUint64(pm_ct) && mxGetNumberOfElements(pm_ct) == 1) {
    uint64_t *ptr = (uint64_t*)mxGetData(pm_ct);
    ct = (CoverTree*)*ptr;
    ct->pm_ctobj = pm;
  }
  mxDestroyArray(pm_ct);
  return ct;
}

static CoverTree* createCoverTree(const mxArray *pm)
{
  CoverTree *ct = (CoverTree*)mxMalloc(sizeof(CoverTree));
  ct->distance =
      getCoverTreeCallback<DistanceCallback>(pm, PROP_DISTANCE_FCN);
  if (ct->distance.function == NULL) {
    mexErrMsgIdAndTxt(COVER_TREE_CLASS":NoDistanceCallback",
                      "Distance callback function must be specified.");
  }
  ct->pre_serialize =
      getCoverTreeCallback<MexCallback>(pm, PROP_PRESERIALIZE_FCN);
  ct->post_deserialize =
      getCoverTreeCallback<MexCallback>(pm, PROP_POSTDESERIALIZE_FCN);
  ct->pm_ctobj = pm;
  return ct;
}

static void makeCoverTreePersistent(CoverTree *ct)
{
  makeCoverTreeCallbackPersistent(ct->distance);
  makeCoverTreeCallbackPersistent(ct->pre_serialize);
  makeCoverTreeCallbackPersistent(ct->post_deserialize);
  mexMakeMemoryPersistent(ct);
}

static mxArray* setCoverTree(const CoverTree *ct)
{
  // Save de CoverTree struct pointer
  mxArray *pm = mxDuplicateArray(ct->pm_ctobj);
  mxArray *pm_cthnd = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  uint64_t *ptr = (uint64_t*)mxGetData(pm_cthnd);
  *ptr = (uint64_t)ct;
  mxSetProperty(pm, 0, PROP_COVER_TREE_HANDLER, pm_cthnd);
  mxDestroyArray(pm_cthnd);
  return pm;
}

static void ctBatchCreate(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
  if (nrhs < 2)
    mexErrMsgTxt("Too few input arguments.");

  if (CT != NULL)
    mexErrMsgTxt("An empty "COVER_TREE_CLASS" object is required.");

  CT = createCoverTree(prhs[0]);
  distance = CT->distance.function;

  v_array<mxArray*> point_set;
  if (! mxIsCell(prhs[1]))
    push(point_set, (mxArray*)prhs[1]);
  else {
    mwSize num_points = mxGetNumberOfElements(prhs[1]);
    alloc(point_set, num_points);
    for(int i = 0; i < num_points; i++)
      push(point_set, mxGetCell(prhs[1], i));
  }

  CT->root_node = batch_create(point_set);
  plhs[0] = setCoverTree(CT);

  mxArray **point_end = point_set.elements + point_set.index;
  for (mxArray **point = point_set.elements; point != point_end; point++)
    mxCreateReference(*point);
  free(point_set);

  persistent::make_persistent<node<mxArray*>* >(persistent::children);
  makeCoverTreePersistent(CT);
  persistent::free(persistent::children);
}

static void insertMxArray(mxArray *pm)
{
  insert(pm, CT->root_node);
  mxCreateReference(pm);
  persistent::make_persistent<node<mxArray*>* >(persistent::children);
}

static void ctInsert(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[])
{
  if (nrhs < 2)
    mexErrMsgTxt("Too few input arguments.");

  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  distance = CT->distance.function;

  if (! mxIsCell(prhs[1]))
    insertMxArray((mxArray*)prhs[1]);
  else {
    mwSize num_points = mxGetNumberOfElements(prhs[1]);
    for(int i = 0; i < num_points; i++)
      insertMxArray(mxGetCell(prhs[1], i));
  }
  persistent::free(persistent::children);
}

static void convertQueryResults(v_array<v_array<d_node<mxArray*> > > vv,
                                mxArray **answers, mxArray **dists)
{
  *answers = mxCreateCellMatrix(vv.index, 1);
  mxArray *temp_dists = mxCreateCellMatrix(vv.index, 1);
  v_array<d_node<mxArray*> > *v_end = vv.elements + vv.index;
  {
    int i = 0;
    for (v_array<d_node<mxArray*> > *v = vv.elements; v != v_end; v++) {
      mxArray *m_ans = mxCreateCellMatrix(v->index, 1);
      mxArray *m_dists = mxCreateDoubleMatrix(v->index, 1, mxREAL);
      double *m_dists_data = mxGetPr(m_dists);
      d_node<mxArray*> *ele_end = v->elements + v->index;
      {
        int j = 0;
        for (d_node<mxArray*> *ele = v->elements; ele != ele_end; ele++) {
          mxSetCell(m_ans, j++, mxCreateReference(ele->n->p));
          *(m_dists_data++) = ele->dist;
        }
      }
      mxSetCell(*answers, i, m_ans);
      mxSetCell(temp_dists, i++, m_dists);
    }
  }
  if (dists != NULL)
    *dists = temp_dists;
  else
    mxDestroyArray(temp_dists);
}

static void ctKNearestNeighbor(int nlhs, mxArray *plhs[],
                               int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3)
    mexErrMsgTxt("Too few input arguments.");

  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  distance = CT->distance.function;

  // Get imput arguments
  CoverTree *query_ct = getCoverTree(prhs[1]);
  if (query_ct == NULL) {
    plhs[0] = mxCreateCellMatrix(0, 0);
    return;
  }
  int k = (int)mxGetScalar(prhs[2]);

  v_array<v_array<d_node<mxArray*> > > results;
  k_nearest_neighbor(CT->root_node, query_ct->root_node, results, k);

  convertQueryResults(results, &plhs[0], (nlhs > 1) ? &plhs[1] : NULL);
  free(results);
}

static void ctEpsilonNearestNeighbor(int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3)
    mexErrMsgTxt("Too few input arguments.");

  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  distance = CT->distance.function;

  // Get imput arguments
  CoverTree *query_ct = getCoverTree(prhs[1]);
  if (query_ct == NULL) {
    plhs[0] = mxCreateCellMatrix(0, 0);
    return;
  }
  double epsilon = mxGetScalar(prhs[2]);

  v_array<v_array<d_node<mxArray*> > > results;
  epsilon_nearest_neighbor(CT->root_node, query_ct->root_node, results,
           epsilon);

  convertQueryResults(results, &plhs[0], (nlhs > 1) ? &plhs[1] : NULL);
  free(results);
}

static void ctUnequalNearestNeighbor(int nlhs, mxArray *plhs[],
                                     int nrhs, const mxArray *prhs[])
{
  if (nrhs < 2)
    mexErrMsgTxt("Too few input arguments.");

  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  distance = CT->distance.function;

  // Get imput arguments
  CoverTree *query_ct = getCoverTree(prhs[1]);
  if (query_ct == NULL) {
    plhs[0] = mxCreateCellMatrix(0, 0);
    return;
  }

  v_array<v_array<d_node<mxArray*> > > results;
  unequal_nearest_neighbor(CT->root_node, query_ct->root_node, results);

  convertQueryResults(results, &plhs[0], (nlhs > 1) ? &plhs[1] : NULL);
  free(results);
}

static void ctDepthDist(int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{
  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  v_array<int> depths;
  depth_dist(CT->root_node.scale, CT->root_node, depths);

  plhs[0] = mxCreateDoubleMatrix(depths.index, 1, mxREAL);
  double *data = mxGetPr(plhs[0]);
  int *depth_end = depths.elements + depths.index;
  for (int *depth = depths.elements; depth != depth_end; depth++, data++) {
    *data = *depth;
  }

  free(depths);
}

static void ctHeightDist(int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[])
{
  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  v_array<int> heights;
  int max_height = height_dist(CT->root_node, heights);

  plhs[0] = mxCreateDoubleMatrix(heights.index, 1, mxREAL);
  double *data = mxGetPr(plhs[0]);
  int *height_end = heights.elements + heights.index;
  for (int *height = heights.elements; height != height_end; height++, data++) {
    *data = *height;
  }

  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleScalar(max_height);
  }

  free(heights);
}

static void ctBreadthDist(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  v_array<int> breadths;
  breadth_dist(CT->root_node, breadths);

  plhs[0] = mxCreateDoubleMatrix(breadths.index, 1, mxREAL);
  double *data = mxGetPr(plhs[0]);
  int *breadth_end = breadths.elements + breadths.index;
  for (int *breadth = breadths.elements; breadth != breadth_end;
        breadth++, data++) {
    *data = *breadth;
  }

  free(breadths);
}

static void deleteChildren(node<mxArray*> *n)
{
  node<mxArray*> *child_end = n->children + n->num_children;
  for (node<mxArray*> *child = n->children; child != child_end; child++) {
    if (child->num_children)
      deleteChildren(child);
    else
      mxDestroyArray(child->p);
  }
  mxFree(n->children);
}

static void ctDelete(int nlhs, mxArray *plhs[],
                     int nrhs, const mxArray *prhs[])
{
  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  deleteChildren(&CT->root_node);
  freeCoverTreeCallback(CT->distance);
  freeCoverTreeCallback(CT->pre_serialize);
  freeCoverTreeCallback(CT->post_deserialize);
  mxFree(CT);
}

template<class Archive>
static void load(Archive &ar, node<mxArray*> &n, unsigned int version)
{
  ar >> n.max_dist;
  ar >> n.parent_dist;
  ar >> n.scale;

  ar >> n.num_children;
  if (n.num_children) {
    n.children =
        (node<mxArray*> *)mxMalloc(sizeof(node<mxArray*>) * n.num_children);
    ar >> boost::serialization::make_array(n.children, n.num_children);
    push(persistent::children, n.children);
    ar >> n.p;
  }
  else {
    n.children = NULL;
    mxArray m;
    ar >> m;
    // Make a copy allocated by MATLAB memory manager
    n.p = mxDuplicateArray(&m);
    ar.reset_object_address(n.p, &m);
    push(persistent::arrays, n.p);
  }
}

template<class Archive>
static void load(Archive &ar, mxArray &p, unsigned int version)
{
  mwSize size;
  ar >> size;

  uint8_t *buffer = (uint8_t*)mxMalloc(size);
  ar >> boost::serialization::make_binary_object(buffer, size);
  mxArray *pm = mxDeserialize(buffer, size);

  if (CT->post_deserialize.function) {
    const mxArray *prhs = (const mxArray*)pm;
    CT->post_deserialize.function(1, &pm, 1, &prhs, CT->post_deserialize.param);
    mxDestroyArray((mxArray*)prhs);
  }

  p = *pm; // pm will be freed when the MEX-function finishes
  mxFree(buffer);
}

static void ctLoad(int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3)
    mexErrMsgTxt("Too few input arguments.");

  if (! mxIsChar(prhs[1]))
    mexErrMsgTxt("The filename argument must be a char array.");
  if (! mxIsChar(prhs[2]))
    mexErrMsgTxt("The mode argument must be a char array.");

  if (CT != NULL)
    mexErrMsgTxt("An empty "COVER_TREE_CLASS" object is required.");

  CT = createCoverTree(prhs[0]);
  char *filename = mxArrayToString(prhs[1]);
  char *mode = mxArrayToString(prhs[2]);

  if (! strcasecmp(mode, "binary")) {
    std::ifstream ifs(filename, std::ios::binary);
    boost::archive::binary_iarchive ar(ifs);
    ar.register_type<mxArray>();
    ar >> CT->root_node;
  }
  else if (! strcasecmp(mode, "text")) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ar(ifs);
    ar.register_type<mxArray>();
    ar >> CT->root_node;
  }
  else {
    mexErrMsgIdAndTxt(COVER_TREE_CLASS":save:UnknownMode", "Unknown mode %s. "
        "Should be 'binary' or 'text'.", mode);
  }
  plhs[0] = setCoverTree(CT);
  mxFree(filename);
  mxFree(mode);

  persistent::make_persistent<mxArray*>(persistent::arrays);
  persistent::make_persistent<node<mxArray*>* >(persistent::children);
  makeCoverTreePersistent(CT);

  persistent::free(persistent::arrays);
  persistent::free(persistent::children);
}

template<class Archive>
static void save(Archive &ar, const node<mxArray*> &n, unsigned int version)
{
  ar << n.max_dist;
  ar << n.parent_dist;
  ar << n.scale;

  ar << n.num_children;
  if (n.num_children) {
    ar << boost::serialization::make_array(n.children, n.num_children);
    ar << n.p;
  }
  else
    ar << *n.p;
}

template<class Archive>
static void save(Archive &ar, const mxArray &p, unsigned int version)
{
  mxArray *psm;

  if (CT->pre_serialize.function) {
    mxArray *plhs = NULL;
    const mxArray *prhs = &p;
    CT->pre_serialize.function(1, &plhs, 1, &prhs, CT->pre_serialize.param);
    psm = mxSerialize(plhs);
    mxDestroyArray(plhs);
  }
  else
    psm = mxSerialize(&p);

  mwSize size = mxGetNumberOfElements(psm);
  ar << size << boost::serialization::make_binary_object(mxGetData(psm), size);
  mxDestroyArray(psm);
}

static void ctSave(int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3)
    mexErrMsgTxt("Too few input arguments.");

  if (! mxIsChar(prhs[1]))
    mexErrMsgTxt("The filename argument must be a char array.");
  if (! mxIsChar(prhs[2]))
    mexErrMsgTxt("The mode argument must be a char array.");

  if (CT == NULL)
    mexErrMsgTxt("A non-empty "COVER_TREE_CLASS" object is required.");

  char *filename = mxArrayToString(prhs[1]);
  char *mode = mxArrayToString(prhs[2]);

  if (! strcasecmp(mode, "binary")) {
    std::ofstream ofs(filename, std::ios::binary);
    boost::archive::binary_oarchive ar(ofs);
    ar.register_type<mxArray>();
    ar << CT->root_node;
  }
  else if (! strcasecmp(mode, "text")) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive ar(ofs);
    ar.register_type<mxArray>();
    ar << CT->root_node;
  }
  else
    mexErrMsgIdAndTxt(COVER_TREE_CLASS":save:UnknownMode", "Unknown mode %s. "
        "Should be 'binary' or 'text'.", mode);

  mxFree(filename);
  mxFree(mode);
}

BOOST_SERIALIZATION_SPLIT_FREE(node<mxArray*>)
BOOST_SERIALIZATION_SPLIT_FREE(mxArray)

// Cleanup function. Invoked by MATLAB before clear the MEX
void cleanup(void)
{
  void **handle_end = PendingHandles + PROP_END;
  for(void **handle = PendingHandles; handle != handle_end; handle++) {
    if (*handle != NULL) {
      dlclose(*handle);
      *handle == NULL;
    }
  }
}

// Entry point to C MEX
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs < 2)
    mexErrMsgTxt("Too few input arguments.");

  if (! mxIsChar(prhs[0]))
    mexErrMsgTxt("Invalid argument #1.");

  // cleanup
  mexAtExit(cleanup);
  cleanup();

  char *str = mxArrayToString(prhs[0]);
  CT = getCoverTree(prhs[1]);

  for (Function *func = functions; func->name != NULL; func++) {
    if (! strcmp(func->name, str)) {
      func->func(nlhs, plhs, nrhs - 1, prhs + 1);
      mxFree(str);
      return;
    }
  }

  mexErrMsgIdAndTxt(COVER_TREE_CLASS":UnknownFunction",
                    "Unknown function name.");
}
