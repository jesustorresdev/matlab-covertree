/* First written by John Langford jl@hunch.net
   Templatization by Dinoj Surendran dinojs@gmail.com
*/

// the files below may not need to be included

/* Whatever structure/class/type is used for P, it must have the following functions defined:

   float distance(P v1, P v2, float upper_bound);
   : this returns the distance between two P objects
   : the distance does not have to be calculated fully if it's more than upper_bound

   v_array<P> parse_points(char *filename);
   : this fills up a v_array of P objects from the input file

   void print(point &P);
   : this prints out the contents of a P object.
*/

#ifndef COVERTREE_H
#define COVERTREE_H

#include <assert.h>
#include <math.h>
#include <values.h>
#include <stdio.h>

#include "stack.h"

// Public interface

template<class P>
struct node {
  P p;
  float max_dist;  // The maximum distance to any grandchild.
  float parent_dist; // The distance to the parent.
  node<P>* children;
  unsigned short int num_children; // The number of children.
  int scale; // Essentially, an upper bound on the distance to any child.
};

template <class P>
struct d_node {
  float dist;
  const node<P> *n;
};

#include "persistent.h"

template<class P>
void print(int depth, node<P> &top_node);

//construction
template<class P>
node<P> new_leaf(const P &p);
template<class P>
node<P> batch_create(v_array<P> points);
template<class P>
void insert(const P &p, node<P> &top_node);
//template<class P>
//void remove(P p, node<P> &top_node); // not yet implemented
//query
template <class P>
void k_nearest_neighbor(const node<P> &top_node, const node<P> &query,
                        v_array<v_array<d_node<P> > > &results, int k);
template <class P>
void epsilon_nearest_neighbor(const node<P> &top_node, const node<P> &query,
                              v_array<v_array<d_node<P> > > &results,
                              float epsilon);
template <class P>
void unequal_nearest_neighbor(const node<P> &top_node, const node<P> &query,
                              v_array<v_array<d_node<P> > > &results);
//information gathering
template <class P>
int height_dist(const node<P> top_node,v_array<int> &heights);
template <class P>
void breadth_dist(const node<P> top_node,v_array<int> &breadths);
template <class P>
void depth_dist(int top_scale, const node<P> top_node,v_array<int> &depths);

//-----------------------------------------------------------------------------

template<class P>
struct ds_node {
  v_array<float> dist;
  P p;
};

float base = 1.3;

float il2 = 1. / log(base);
inline float dist_of_scale (int s)
{
  return pow(base, s);
}

inline int get_scale(float d)
{
  return (int) ceilf(il2 * log(d));
}

int min(int f1, int f2)
{
  if ( f1 <= f2 )
    return f1;
  else
    return f2;
}

template<class P>
P max(P f1, P f2)
{
  if ( f1 <= f2 )
    return f2;
  else
    return f1;
}

template<class P>
node<P> new_node(const P &p)
{
  node<P> new_node;
  new_node.p = p;
  return new_node;
}

template<class P>
node<P> new_leaf(const P &p)
{
  node<P> new_leaf = {p,0.,0.,NULL,0,INT_MIN};
  return new_leaf;
}

template<class P>
d_node<P> new_dnode(const node<P> *n, float dist = 0.)
{
  d_node<P> new_dnode = {dist, n};
  return new_dnode;
}

template<class P>
float max_set(v_array<ds_node<P> > &v)
{
  float max = 0.;
  for (int i = 0; i < v.index; i++)
    if ( max < v[i].dist.last())
      max = v[i].dist.last();
  return max;
}

void print_space(int s)
{
  for (int i = 0; i < s; i++)
    printf(" ");
}

template<class P>
void print(int depth, node<P> &top_node)
{
  print_space(depth);
  print(top_node.p);
  if ( top_node.num_children > 0 ) {
    print_space(depth); printf("scale = %i\n",top_node.scale);
    print_space(depth); printf("max_dist = %f\n",top_node.max_dist);
    print_space(depth); printf("num children = %i\n",top_node.num_children);
    for (int i = 0; i < top_node.num_children;i++)
      print(depth+1, top_node.children[i]);
  }
}

template<class P>
void split(v_array<ds_node<P> >& point_set,
           v_array<ds_node<P> >& far_set,
           int max_scale)
{
  unsigned int new_index = 0;
  float fmax = dist_of_scale(max_scale);
  for (int i = 0; i < point_set.index; i++){
    if (point_set[i].dist.last() <= fmax) {
      point_set[new_index++] = point_set[i];
    }
    else
      push(far_set,point_set[i]);
  }
  point_set.index=new_index;
}

template<class P>
void dist_split(v_array<ds_node<P> >& point_set,
                v_array<ds_node<P> >& new_point_set,
                P new_point,
                int max_scale)
{
  unsigned int new_index = 0;
  float fmax = dist_of_scale(max_scale);
  for(int i = 0; i < point_set.index; i++)
    {
      float new_d;
      new_d = distance(new_point, point_set[i].p, fmax);
      if (new_d <= fmax ) {
        push(point_set[i].dist, new_d);
        push(new_point_set,point_set[i]);
      }
      else
        point_set[new_index++] = point_set[i];
    }
  point_set.index = new_index;
}

/*
   max_scale is the maximum scale of the node we might create here.
   point_set contains points which are 2*max_scale or less away.
*/

template <class P>
node<P> batch_insert(const P& p,
                     int max_scale,
                     v_array<ds_node<P> >& point_set,
                     v_array<ds_node<P> >& consumed_set,
                     v_array<v_array<ds_node<P> > >& stack)
{
  if (point_set.index == 0)
    return new_leaf(p);
  else {
    float max_dist = max_set(point_set); //O(|point_set|)
    int next_scale = min (max_scale - 1, get_scale(max_dist));
    if (next_scale == INT_MIN) // We have points with distance 0.
      {
        v_array<node<P> > children;
        push(children,new_leaf(p));
        while (point_set.index > 0)
          {
            push(children,new_leaf(point_set.last().p));
            push(consumed_set,point_set.last());
            point_set.decr();
          }
        node<P> n = new_node(p);
        n.scale = INT_MIN; // A magic number meant to be larger than all scales.
        n.max_dist = 0;
        alloc(children,children.index);
        n.num_children = children.index;
        n.children = children.elements;
        push(persistent::children, n.children);
        return n;
      }
    else
      {
        v_array<ds_node<P> > far = pop(stack);
        split(point_set,far,max_scale); //O(|point_set|)

        node<P> child = batch_insert(p, next_scale, point_set, consumed_set, stack);

        if (point_set.index == 0)
          {
            push(stack,point_set);
            point_set=far;
            return child;
          }
        else {
          node<P> n = new_node(p);
          v_array<node<P> > children;
          push(children, child);
          v_array<ds_node<P> > new_point_set = pop(stack);
          v_array<ds_node<P> > new_consumed_set = pop(stack);
          while (point_set.index != 0) { //O(|point_set| * num_children)
            P new_point = point_set.last().p;
            float new_dist = point_set.last().dist.last();
            push(consumed_set, point_set.last());
            point_set.decr();

            dist_split(point_set, new_point_set, new_point, max_scale); //O(|point_saet|)
            dist_split(far,new_point_set,new_point,max_scale); //O(|far|)

            node<P> new_child =
              batch_insert(new_point, next_scale, new_point_set, new_consumed_set, stack);
            new_child.parent_dist = new_dist;

            push(children, new_child);

            float fmax = dist_of_scale(max_scale);
            for(int i = 0; i< new_point_set.index; i++) //O(|new_point_set|)
              {
                new_point_set[i].dist.decr();
                if (new_point_set[i].dist.last() <= fmax)
                  push(point_set, new_point_set[i]);
                else
                  push(far, new_point_set[i]);
              }
            for(int i = 0; i< new_consumed_set.index; i++) //O(|new_point_set|)
              {
                new_consumed_set[i].dist.decr();
                push(consumed_set, new_consumed_set[i]);
              }
            new_point_set.index = 0;
            new_consumed_set.index = 0;
          }
          push(stack,new_point_set);
          push(stack,new_consumed_set);
          push(stack,point_set);
          point_set=far;
          n.scale = max_scale;
          n.max_dist = max_set(consumed_set);
          alloc(children,children.index);
          n.num_children = children.index;
          n.children = children.elements;
          push(persistent::children, n.children);
          return n;
        }
      }
  }
}

template<class P>
node<P> batch_create(v_array<P> points)
{
  assert(points.index > 0);
  v_array<ds_node<P> > point_set;
  v_array<v_array<ds_node<P> > > stack;

  for (int i = 1; i < points.index; i++) {
    ds_node<P> temp;
    push(temp.dist, distance(points[0], points[i], MAXFLOAT));
    temp.p = points[i];
    push(point_set,temp);
  }

  v_array<ds_node<P> > consumed_set;

  float max_dist = max_set(point_set);

  node<P> top = batch_insert (points[0],
                           get_scale(max_dist),
                            point_set,
                            consumed_set,
                            stack);
  for (int i = 0; i<consumed_set.index;i++)
    free(consumed_set[i].dist);
  free(consumed_set);
  for (int i = 0; i<stack.index;i++)
    free(stack[i]);
  free(stack);
  free(point_set);
  return top;
}

void add_height(int d, v_array<int> &heights)
{
  if (heights.index <= d)
    for(;heights.index <= d;)
      push(heights,0);
  heights[d] = heights[d] + 1;
}

template <class P>
int height_dist(const node<P> top_node, v_array<int> &heights)
{
  if (top_node.num_children == 0)
    {
      add_height(0,heights);
      return 0;
    }
  else
    {
      int max_v=0;
      for (int i = 0; i<top_node.num_children ;i++)
        {
          int d = height_dist(top_node.children[i], heights);
          if (d > max_v)
            max_v = d;
        }
      add_height(1 + max_v, heights);
      return (1 + max_v);
    }
}

template <class P>
void depth_dist(int top_scale, const node<P> top_node,v_array<int> &depths)
{
  if (top_node.num_children > 0)
      for (int i = 0; i<top_node.num_children ;i++)
        {
          add_height(top_scale - top_node.scale, depths);
          depth_dist(top_scale, top_node.children[i], depths);
        }
}

template <class P>
void breadth_dist(const node<P> top_node, v_array<int> &breadths)
{
  if (top_node.num_children == 0)
    add_height(0,breadths);
  else
    {
      for (int i = 0; i<top_node.num_children ;i++)
        breadth_dist(top_node.children[i], breadths);
      add_height(top_node.num_children, breadths);
    }
}

template <class P>
inline float compare(const d_node<P> *p1, const d_node<P>* p2)
{
  return p1 -> dist - p2 -> dist;
}

template <class P>
inline void SWAP (d_node<P>* a, d_node<P>* b)
{
      d_node<P> tmp = *a;
      *a = *b;
      *b = tmp;
}

template <class P>
void halfsort (v_array<d_node<P> > cover_set)
{

  if (cover_set.index <= 1)
    return;
  register d_node<P> *base_ptr =  cover_set.elements;

  d_node<P> *hi = &base_ptr[cover_set.index - 1];
  d_node<P> *right_ptr = hi;
  d_node<P> *left_ptr;

  while (right_ptr > base_ptr)
    {
      d_node<P> *mid = base_ptr + ((hi - base_ptr) >> 1);

      if (compare ( mid,  base_ptr) < 0.)
        SWAP (mid, base_ptr);
      if (compare ( hi,  mid) < 0.)
        SWAP (mid, hi);
      else
        goto jump_over;
      if (compare ( mid,  base_ptr) < 0.)
        SWAP (mid, base_ptr);
    jump_over:;

      left_ptr  = base_ptr + 1;
      right_ptr = hi - 1;

      do
        {
          while (compare (left_ptr, mid) < 0.)
            left_ptr ++;

          while (compare (mid, right_ptr) < 0.)
            right_ptr --;

          if (left_ptr < right_ptr)
            {
              SWAP (left_ptr, right_ptr);
              if (mid == left_ptr)
                mid = right_ptr;
              else if (mid == right_ptr)
                mid = left_ptr;
              left_ptr ++;
              right_ptr --;
            }
          else if (left_ptr == right_ptr)
            {
              left_ptr ++;
              right_ptr --;
              break;
            }
        }
      while (left_ptr <= right_ptr);

      hi = right_ptr;
    }
}

template <class P>
v_array<v_array<d_node<P> > >
    get_cover_sets(v_array<v_array<v_array<d_node<P> > > > &spare_cover_sets)
{
  v_array<v_array<d_node<P> > > ret = pop(spare_cover_sets);
  while (ret.index < 101)
    {
      v_array<d_node<P> > temp;
      push(ret, temp);
    }
  return ret;
}

inline bool shell(float parent_query_dist, float child_parent_dist,
                  float upper_bound)
{
  return parent_query_dist - child_parent_dist <= upper_bound;
  //    && child_parent_dist - parent_query_dist <= upper_bound;
}

int internal_k =1;
void update_k(float *k_upper_bound, float upper_bound)
{
  float *end = k_upper_bound + internal_k-1;
  float *begin = k_upper_bound;
  for (;end != begin; begin++)
    {
      if (upper_bound < *(begin+1))
        *begin = *(begin+1);
      else {
        *begin = upper_bound;
        break;
      }
    }
  if (end == begin)
    *begin = upper_bound;
}
float *alloc_k()
{
  return (float *)mxMalloc(sizeof(float) * internal_k);
}
void set_k(float* begin, float max)
{
  for(float *end = begin+internal_k;end != begin; begin++)
    *begin = max;
}

float internal_epsilon =0.;
void update_epsilon(float *upper_bound, float new_dist) {}
float *alloc_epsilon()
{
  return (float *)mxMalloc(sizeof(float));
}
void set_epsilon(float* begin, float max)
{
  *begin = internal_epsilon;
}

void update_unequal(float *upper_bound, float new_dist)
{
  if (new_dist != 0.)
    *upper_bound = new_dist;
}
float* (*alloc_unequal)() = alloc_epsilon;
void set_unequal(float* begin, float max)
{
  *begin = max;
}

void (*update)(float *foo, float bar) = update_k;
void (*setter)(float *foo, float bar) = set_k;
float* (*alloc_upper)() = alloc_k;

template <class P>
inline void copy_zero_set(node<P>* query_chi,
                          float* new_upper_bound,
                          v_array<d_node<P> > &zero_set,
                          v_array<d_node<P> > &new_zero_set)
{
  new_zero_set.index = 0;
  d_node<P> *end = zero_set.elements + zero_set.index;
  for (d_node<P> *ele = zero_set.elements; ele != end ; ele++)
    {
      float upper_dist = *new_upper_bound + query_chi->max_dist;
      if (shell(ele->dist, query_chi->parent_dist, upper_dist))
        {
          float d = distance(query_chi->p, ele->n->p, upper_dist);

          if (d <= upper_dist)
            {
              if (d < *new_upper_bound)
                update(new_upper_bound, d);
              d_node<P> temp = {d, ele->n};
              push(new_zero_set,temp);
            }
        }
    }
}

template <class P>
inline void copy_cover_sets(node<P>* query_chi,
                            float* new_upper_bound,
                            v_array<v_array<d_node<P> > > &cover_sets,
                            v_array<v_array<d_node<P> > > &new_cover_sets,
                            int current_scale,
                            int max_scale)
{
  for (; current_scale <= max_scale; current_scale++)
    {
      d_node<P>* ele = cover_sets[current_scale].elements;
      d_node<P>* end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
      for (; ele != end; ele++)
        {
          float upper_dist = *new_upper_bound + query_chi->max_dist + ele->n->max_dist;
          if (shell(ele->dist, query_chi->parent_dist, upper_dist))
            {
              float d = distance(query_chi->p, ele->n->p, upper_dist);

              if (d <= upper_dist)
                {
                  if (d < *new_upper_bound)
                    update(new_upper_bound,d);
                  d_node<P> temp = {d, ele->n};
                  push(new_cover_sets[current_scale],temp);
                }
            }
        }
    }
}

template <class P>
void print_query(const node<P> *top_node)
{
  printf ("query = \n");
  print(top_node->p);
  if ( top_node->num_children > 0 ) {
    printf("scale = %i\n",top_node->scale);
    printf("max_dist = %f\n",top_node->max_dist);
    printf("num children = %i\n",top_node->num_children);
  }
}

template <class P>
void print_cover_sets(v_array<v_array<d_node<P> > > &cover_sets,
                      v_array<d_node<P> > &zero_set,
                      int current_scale,
                      int max_scale)
{
  printf("cover set = \n");
  for (; current_scale <= max_scale; current_scale++)
    {
      d_node<P> *ele = cover_sets[current_scale].elements;
      d_node<P> *end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
      printf ("%i\n", current_scale);
      for (; ele != end; ele++)
        {
          node<P> *n = (node<P> *)ele->n;
          print(n->p);
        }
    }
  d_node<P> *end = zero_set.elements + zero_set.index;
  printf ("infinity\n");
  for (d_node<P> *ele = zero_set.elements; ele != end ; ele++)
    {
      node<P> *n = (node<P> *)ele->n;
      print(n->p);
    }
}

/*
  An optimization to consider:
  Make all distance evaluations occur in descend.

  Instead of passing a cover_set, pass a stack of cover sets.  The
  last element holds d_nodes with your distance.  The next lower
  element holds a d_node with the distance to your query parent,
  next = query grand parent, etc..

  Compute distances in the presence of the tighter upper bound.
 */

template <class P>
inline
void descend(const node<P>* query,
             float* upper_bound,
             int current_scale,
             int &max_scale,
             int top_scale,
             v_array<v_array<d_node<P> > > &cover_sets,
             v_array<d_node<P> > &zero_set)
{
  d_node<P> *end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
  for (d_node<P> *parent = cover_sets[current_scale].elements; parent != end; parent++)
    {
      const node<P> *par = parent->n;
      float upper_dist = *upper_bound + query->max_dist + query->max_dist;
      if (parent->dist <= upper_dist + par->max_dist)
        {
          node<P> *chi = par->children;
          if (parent->dist <= upper_dist + chi->max_dist)
            if (chi->num_children > 0)
              {
                int chi_scale = top_scale - chi->scale;
                if (max_scale < chi_scale)
                  max_scale = chi_scale;
                d_node<P> temp = {parent->dist, chi};
                push(cover_sets[chi_scale], temp);
              }
            else if (parent->dist <= upper_dist)
              {
                d_node<P> temp = {parent->dist, chi};
                push(zero_set, temp);
              }
          node<P> *child_end = par->children + par->num_children;
          for (chi++; chi != child_end; chi++)
            {
              float upper_chi = *upper_bound + chi->max_dist + query->max_dist + query->max_dist;
              if (shell(parent->dist, chi->parent_dist, upper_chi))
                {
                  float d = distance(query->p, chi->p, upper_chi);
                  if (d <= upper_chi)
                    {
                      if (d < *upper_bound)
                        update(upper_bound, d);
                      if (chi->num_children > 0)
                        {
                          int chi_scale = top_scale - chi->scale;
                          if (max_scale < chi_scale)
                            max_scale = chi_scale;
                          d_node<P> temp = {d, chi};
                          push(cover_sets[chi_scale],temp);
                        }
                      else
                        if (d <= upper_chi - chi->max_dist)
                          {
                            d_node<P> temp = {d, chi};
                            push(zero_set, temp);
                          }
                    }
                }
            }
        }
    }
}

template <class P>
void brute_nearest(const node<P>* query,v_array<d_node<P> > zero_set,
                   float* upper_bound,
                   v_array<v_array<d_node <P> > > &results,
                   v_array<v_array<d_node<P> > > &spare_zero_sets)
{
  if (query->num_children > 0)
    {
      v_array<d_node<P> > new_zero_set = pop(spare_zero_sets);
      node<P> * query_chi = query->children;
      brute_nearest(query_chi, zero_set, upper_bound, results, spare_zero_sets);
      float* new_upper_bound = alloc_upper();

      node<P> *child_end = query->children + query->num_children;
      for (query_chi++;query_chi != child_end; query_chi++)
        {
          setter(new_upper_bound,*upper_bound + query_chi->parent_dist);
          copy_zero_set(query_chi, new_upper_bound, zero_set, new_zero_set);
          brute_nearest(query_chi, new_zero_set, new_upper_bound, results, spare_zero_sets);
        }
      mxFree(new_upper_bound);
      new_zero_set.index = 0;
      push(spare_zero_sets, new_zero_set);
    }
  else
    {
      v_array<d_node<P> > temp;
      push(temp, new_dnode(query));
      d_node<P> *end = zero_set.elements + zero_set.index;
      for (d_node<P> *ele = zero_set.elements; ele != end ; ele++)
        if (ele->dist <= *upper_bound)
            push(temp, new_dnode(ele->n, ele->dist));
      push(results, temp);
    }
}

template <class P>
void internal_batch_nearest_neighbor(const node<P> *query,
                                     v_array<v_array<d_node<P> > > &cover_sets,
                                     v_array<d_node<P> > &zero_set,
                                     int current_scale,
                                     int max_scale,
                                     int top_scale,
                                     float* upper_bound,
                                     v_array<v_array<d_node<P> > > &results,
                                     v_array<v_array<v_array<d_node<P> > > > &spare_cover_sets,
                                     v_array<v_array<d_node<P> > > &spare_zero_sets)
{
  if (current_scale > max_scale) // All remaining points are in the zero set.
    brute_nearest(query, zero_set, upper_bound, results, spare_zero_sets);
  else
    if (query->scale >= (top_scale - current_scale) && query->scale != INT_MIN)
      // Our query has too much scale.  Reduce.
      {
        node<P> *query_chi = query->children;
        v_array<d_node<P> > new_zero_set = pop(spare_zero_sets);
        v_array<v_array<d_node<P> > > new_cover_sets = get_cover_sets(spare_cover_sets);
        float* new_upper_bound = alloc_upper();

        node<P> *child_end = query->children + query->num_children;
        for (query_chi++; query_chi != child_end; query_chi++)
          {
            setter(new_upper_bound,*upper_bound + query_chi->parent_dist);
            copy_zero_set(query_chi, new_upper_bound, zero_set, new_zero_set);
            copy_cover_sets(query_chi, new_upper_bound, cover_sets, new_cover_sets,
                              current_scale, max_scale);
            internal_batch_nearest_neighbor(query_chi, new_cover_sets, new_zero_set,
                                            current_scale, max_scale, top_scale,
                                            new_upper_bound, results, spare_cover_sets,
                                            spare_zero_sets);
          }
        mxFree(new_upper_bound);
        new_zero_set.index = 0;
        push(spare_zero_sets, new_zero_set);
        push(spare_cover_sets, new_cover_sets);
        internal_batch_nearest_neighbor(query->children, cover_sets, zero_set,
                                        current_scale, max_scale, top_scale,
                                        upper_bound, results, spare_cover_sets,
                                        spare_zero_sets);
      }
    else // reduce cover set scale
      {
        halfsort(cover_sets[current_scale]);
        descend(query, upper_bound, current_scale, max_scale, top_scale,
                cover_sets, zero_set);
        cover_sets[current_scale++].index = 0;
        internal_batch_nearest_neighbor(query, cover_sets, zero_set,
                                        current_scale, max_scale, top_scale,
                                        upper_bound, results, spare_cover_sets,
                                        spare_zero_sets);
      }
}

template <class P>
void batch_nearest_neighbor(const node<P> &top_node,
                            const node<P> &query,
                            v_array<v_array<d_node<P> > > &results)
{
  v_array<v_array<v_array<d_node<P> > > > spare_cover_sets;
  v_array<v_array<d_node<P> > > spare_zero_sets;

  v_array<v_array<d_node<P> > > cover_sets = get_cover_sets(spare_cover_sets);
  v_array<d_node<P> > zero_set = pop(spare_zero_sets);

  float* upper_bound = alloc_upper();
  setter(upper_bound,MAXFLOAT);

  float top_dist = distance(query.p, top_node.p, MAXFLOAT);
  update(upper_bound, top_dist);

  d_node<P> temp = {top_dist, &top_node};
  push(cover_sets[0], temp);

  internal_batch_nearest_neighbor(&query, cover_sets, zero_set, 0, 0,
                                  top_node.scale, upper_bound, results,
                                  spare_cover_sets, spare_zero_sets);

  mxFree(upper_bound);
  push(spare_cover_sets, cover_sets);

  for (int i = 0; i < spare_cover_sets.index; i++)
    {
      v_array<v_array<d_node<P> > > cover_sets = spare_cover_sets[i];
      for (int j = 0; j < cover_sets.index; j++)
        free(cover_sets[j]);
      free(cover_sets);
    }
  free(spare_cover_sets);

  push(spare_zero_sets, zero_set);

  for (int i = 0; i < spare_zero_sets.index; i++)
    free(spare_zero_sets[i]);
  free(spare_zero_sets);
}

template <class P>
void k_nearest_neighbor(const node<P> &top_node,
                        const node<P> &query,
                        v_array<v_array<d_node<P> > > &results,
                        int k)
{
  internal_k = k;
  update = update_k;
  setter = set_k;
  alloc_upper = alloc_k;

  batch_nearest_neighbor(top_node, query,results);
}

template <class P>
void epsilon_nearest_neighbor(const node<P> &top_node,
                              const node<P> &query,
                              v_array<v_array<d_node<P> > > &results,
                              float epsilon)
{
  internal_epsilon = epsilon;
  update = update_epsilon;
  setter = set_epsilon;
  alloc_upper = alloc_epsilon;

  batch_nearest_neighbor(top_node, query,results);
}

template <class P>
void unequal_nearest_neighbor(const node<P> &top_node,
                              const node<P> &query,
                              v_array<v_array<d_node<P> > > &results)
{
  update = update_unequal;
  setter = set_unequal;
  alloc_upper = alloc_unequal;

  batch_nearest_neighbor(top_node, query, results);
}

template <class P>
void add_children(node<P> &par, const node<P> &n)
{
  par.children = (node<P> *)mxRealloc(par.children,
                                      sizeof(node<P>) * (++par.num_children));
  par.children[par.num_children-1] = n;
}

template <class P>
inline
void descend_insert(const P &p,
                    float* upper_bound,
                    int current_scale,
                    int &max_scale,
                    int top_scale,
                    v_array<v_array<d_node<P> > > &cover_sets,
                    v_array<d_node<P> > &parent_set)
{
  *upper_bound = dist_of_scale(top_scale - current_scale);
  d_node<P> *end = cover_sets[current_scale].elements + cover_sets[current_scale].index;
  for (d_node<P> *parent = cover_sets[current_scale].elements; parent != end; parent++)
    {
      const node<P> *par = parent->n;
      node<P> *chi = par->children;
      if (parent->dist <= *upper_bound)
        if (chi->scale != INT_MIN && parent->dist <= chi->max_dist)
          {
            int chi_scale = top_scale - chi->scale;
            if (max_scale < chi_scale)
              max_scale = chi_scale;
            d_node<P> temp = {parent->dist, chi};
            push(cover_sets[chi_scale], temp);
          }
        else
          {
            d_node<P> temp = {parent->dist, chi};
            push(parent_set, temp);
          }
      node<P> *child_end = par->children + par->num_children;
      for (chi++; chi != child_end; chi++)
        {
          if (shell(parent->dist, chi->parent_dist, *upper_bound))
            {
              float d = distance(p, chi->p, *upper_bound);
              if (d <= *upper_bound)
                if (chi->scale != INT_MIN && d <= chi->max_dist)
                  {
                    int chi_scale = top_scale - chi->scale;
                    if (max_scale < chi_scale)
                      max_scale = chi_scale;
                    d_node<P> temp = {d, chi};
                    push(cover_sets[chi_scale], temp);
                  }
                else
                  {
                    d_node<P> temp = {d, chi};
                    push(parent_set, temp);
                  }
            }
        }
    }
}

template <class P>
void internal_insert(const P &p,
                     v_array<v_array<d_node<P> > > &cover_sets,
                     v_array<d_node<P> > &parent_set,
                     int current_scale,
                     int max_scale,
                     int top_scale,
                     float* upper_bound)
{
  if (current_scale > max_scale)
    {
      d_node<P> *end = parent_set.elements + parent_set.index;
      for (d_node<P> *ele = parent_set.elements; ele != end; ele++)
        if (ele->dist <= *upper_bound)
          {
            node<P> n = new_leaf(p);
            n.parent_dist = ele->dist;
            node<P> *par = (node<P>*)ele->n;
            if (par->num_children != 0 && ele->dist == 0.)
              add_children(*par, n);
            else
              {
                v_array<node<P> > children;
                alloc(children, 2);
                push(children, *par);
                children[0].parent_dist = 0.;
                push(children, n);
                par->max_dist = ele->dist;
                par->scale = get_scale(ele->dist);
                par->num_children = children.index;
                par->children = children.elements;
              }
            push(persistent::children, par->children);
            return;
          }

      d_node<P> *parent = cover_sets[max_scale].elements;
      node<P> *par = (node<P>*)parent->n;
      node<P> n = new_leaf(p);
      n.parent_dist = parent->dist;
      add_children(*par, n);
      push(persistent::children, par->children);
    }
  else
    {
      descend_insert(p, upper_bound, current_scale, max_scale, top_scale,
                     cover_sets, parent_set);
      internal_insert(p, cover_sets, parent_set, ++current_scale, max_scale,
                      top_scale, upper_bound);
    }
}

template <class P>
void insert(const P &p, node<P> &top_node)
{
  float top_dist = distance(p, top_node.p, MAXFLOAT);
  if (top_node.scale != 100 && top_dist <= top_node.max_dist)
    {
      v_array<v_array<v_array<d_node<P> > > > spare_cover_sets;
      v_array<v_array<d_node<P> > > cover_sets = get_cover_sets(spare_cover_sets);
      v_array<d_node<P> > parent_set;

      d_node<P> temp = {top_dist, &top_node};
      push(cover_sets[0], temp);

      float upper_bound = top_dist;
      internal_insert(p, cover_sets, parent_set, 0, 0, top_node.scale,
                      &upper_bound);

      free(cover_sets);
      free(parent_set);
    }
  else if (top_node.num_children != 0 && top_dist == 0.)
    {
      node<P> n = new_leaf(p);
      n.parent_dist = 0.;
      add_children(top_node, n);
      push(persistent::children, top_node.children);
    }
  else
    {
      v_array<node<P> > children;
      alloc(children, 2);
      push(children, top_node);
      children[0].parent_dist = 0.;
      node<P> n = new_leaf(p);
      n.parent_dist = top_dist;
      push(children, n);
      top_node.max_dist = top_dist;
      top_node.scale = get_scale(top_dist);
      top_node.num_children = children.index;
      top_node.children = children.elements;
      push(persistent::children, top_node.children);
    }
}

#endif
