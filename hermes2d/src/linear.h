// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_LINEAR_H
#define __H2D_LINEAR_H

#include "common.h"
#include "solution.h"


const double H2D_EPS_LOW    = 0.0014;
const double H2D_EPS_NORMAL = 0.0008;
const double H2D_EPS_HIGH   = 0.0003;


/// Linearizer is a utility class which converts a higher-order FEM solution defined on
/// a curvilinear, irregular mesh to a linear FEM solution defined on a straight-edged,
/// regular mesh. This is done by adaptive refinement of the higher-order mesh and its
/// subsequent regularization. The linearized mesh can then be easily displayed or
/// exported to standard formats. The class correctly handles discontinuities in the
/// solution (e.g., gradients or in Hcurl) by inserting double vertices where necessary.
/// Linearizer also serves as a container for the resulting linearized mesh.
///
class H2D_API Linearizer // (implemented in linear1.cpp)
{
public:

  Linearizer();
  ~Linearizer();

  void process_solution(MeshFunction* sln, int item = H2D_FN_VAL_0,
                        double eps = H2D_EPS_NORMAL, double max_abs = -1.0,
                        MeshFunction* xdisp = NULL, MeshFunction* ydisp = NULL,
                        double dmult = 1.0);

  void lock_data() const { pthread_mutex_lock(&data_mutex); }
  void unlock_data() const { pthread_mutex_unlock(&data_mutex); }

  double3* get_vertices() const { return verts; }
  int get_num_vertices() const { return nv; }

  int3* get_triangles() const { return tris; }
  int get_num_triangles() const { return nt; }

  int3* get_edges() const { return edges; }
  int get_num_edges() const { return ne; }

  double get_min_value() const { return min_val; }
  double get_max_value() const { return max_val; }
  virtual void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

  virtual void save_data(const char* filename);
  virtual void load_data(const char* filename);

  void free();

protected:

  MeshFunction* sln;
  int item, ia, ib;

  double eps, max, cmax;
  bool auto_max;

  MeshFunction *xdisp, *ydisp;
  double dmult;

  double3* verts;  ///< vertices: (x, y, value) triplets
  int4* info;      ///< info[0]=p1, info[1]=p2, info[2]=next vertex in hash
  int3* tris;      ///< triangles: vertex index triplets
  int3* edges;     ///< edges: pairs of vertex indices
  int* hash_table; ///< hash table

  int nv, nt, ne; ///< numbers of vertices, triangles and edges
  int cv, ct, ce; ///< capacities of vertex, triangle and edge arrays
  int del_slot;   ///< free slot index after a triangle which was deleted
  int mask;       ///< hash table mask = size-1

  bool curved, disp;
  double min_val, max_val;

  int get_vertex(int p1, int p2, double x, double y, double value);
  int get_top_vertex(int id, double value);
  int peek_vertex(int p1, int p2);

  int hash(int p1, int p2) { return (984120265*p1 + 125965121*p2) & mask; }

  int add_vertex()
  {
    if (nv >= cv)
    {
      cv *= 2;
      verts = (double3*) realloc(verts, sizeof(double3) * cv);
      info = (int4*) realloc(info, sizeof(int4) * cv);
      verbose("Linearizer::add_vertex(): realloc to %d", cv);
    }
    return nv++;
  }

  void add_triangle(int iv0, int iv1, int iv2);

  void del_triangle(int index)
  {
    del_slot = index;
  }

  void add_edge(int iv1, int iv2, int marker)
  {
    if (ne >= ce) edges = (int3*) realloc(edges, sizeof(int3) * (ce = ce * 3 / 2));
    edges[ne][0] = iv1;
    edges[ne][1] = iv2;
    edges[ne++][2] = marker;
  }

  void get_gv_a_b(int mask, int& a, int& b)
  {
    a = b = 0;
    if (mask >= 0x40) { a = 1; mask >>= 6; }
    while (!(mask & 1)) { mask >>= 1; b++; }
  }

  void process_triangle(int iv0, int iv1, int iv2, int level,
                        scalar* val, double* phx, double* phy, int* indices);

  void process_quad(int iv0, int iv1, int iv2, int iv3, int level,
                    scalar* val, double* phx, double* phy, int* indices);

  void process_edge(int iv1, int iv2, int marker);
  void regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2);
  void find_min_max();
  void print_hash_stats();

  mutable pthread_mutex_t data_mutex;

  static void calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y); ///< Calculates AABB from an array of X-axis and Y-axis coordinates. The distance between values in the array is stride bytes.
};


/// Like the Linearizer, but generates a triangular mesh showing polynomial
/// orders in a space, hence the funky name.
///
class H2D_API Orderizer : public Linearizer // (implemented in linear2.cpp)
{
public:

  Orderizer();
  ~Orderizer();

  void process_solution(Space* space);

  int get_labels(int*& lvert, char**& ltext, double2*& lbox) const
        { lvert = this->lvert; ltext = this->ltext; lbox = this->lbox; return nl; };

  virtual void save_data(const char* filename);
  virtual void load_data(const char* filename);

protected:

  char  buffer[1000];
  char* labels[11][11];

  int  nl, cl1, cl2, cl3;
  int* lvert;
  char** ltext;
  double2* lbox;

};


/// "Vectorizer" is a Linearizer for vector solutions. The only difference is
/// that linearized vertices are vector-valued. Also, regularization of the
/// resulting mesh is not attempted. The class can handle different meshes in
/// both X and Y components.
///
class H2D_API Vectorizer : public Linearizer // (implemented in linear3.cpp)
{
public:

  Vectorizer();
  ~Vectorizer();

  void process_solution(MeshFunction* xsln, int xitem, MeshFunction* ysln, int yitem, double eps);

public: //accessors
  double4* get_vertices() const { return verts; }
  int get_num_vertices() const { return nv; }

  int2* get_dashes() const { return dashes; }
  int get_num_dashes() const { return nd; }

  double get_min_value() const { return min_val; }
  double get_max_value() const { return max_val; }

  virtual void save_data(const char* filename);
  virtual void load_data(const char* filename);
  virtual void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

  void free();

protected:

  MeshFunction *xsln, *ysln;
  int xitem, yitem;
  int xia, xib, yia, yib;
  double4* verts;  ///< vertices: (x, y, xvalue, yvalue) quadruples
  int2* dashes;
  int nd, ed, cd;

  int get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue);
  int create_vertex(double x, double y, double xvalue, double yvalue);
  void process_dash(int iv1, int iv2);

  int add_vertex()
  {
    if (nv >= cv)
    {
      cv *= 2;
      verts = (double4*) realloc(verts, sizeof(double4) * cv);
      info = (int4*) realloc(info, sizeof(int4) * cv);
      verbose("Vectorizer::add_vertex(): realloc to %d", cv);
    }
    return nv++;
  }

  void add_dash(int iv1, int iv2)
  {
    if (nd >= cd) dashes = (int2*) realloc(dashes, sizeof(int2) * (cd = cd * 3 / 2));
    dashes[nd][0] = iv1;
    dashes[nd++][1] = iv2;
  }

  void push_transform(int son)
  {
    xsln->push_transform(son);
    if (ysln != xsln) ysln->push_transform(son);
  }

  void pop_transform()
  {
    xsln->pop_transform();
    if (ysln != xsln) ysln->pop_transform();
  }

  void process_triangle(int iv0, int iv1, int iv2, int level,
                        scalar* xval, scalar* yval, double* phx, double* phy, int* indices);

  void process_quad(int iv0, int iv1, int iv2, int iv3, int level,
                    scalar* xval, scalar* yval, double* phx, double* phy, int* indices);

  void find_min_max();

};


// maximum subdivision level (2^N)
const int LIN_MAX_LEVEL = 6;


#define lin_init_array(array, type, c, e) \
  if (c < e) { \
    if (array != NULL) ::free(array); \
    array = (type*) malloc(sizeof(type) * (c = e)); }

#define lin_free_array(array, n, c) \
  if (array != NULL) { \
    ::free(array); array = NULL; \
    n = c = 0; }



#endif
