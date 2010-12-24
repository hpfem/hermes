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

#include "h2d_common.h"
#include "linear.h"
#include "refmap.h"


//// linearization "quadrature" ////////////////////////////////////////////////////////////////////

// The tables with index zero are for obtaining solution values at the element
// vertices. Index one tables serve for the retrieval of interior values. Index one tables
// are used for adaptive approximation of the solution by transforming their points to sub-elements.
// Actually, the tables contain two levels of refinement -- this is an optimization to reduce
// the number of calls to sln->get_values().

static double3 lin_pts_0_tri[] =
{
  { -1.0, -1.0, 0.0 },
  {  1.0, -1.0, 0.0 },
  { -1.0,  1.0, 0.0 }
};

static double3 lin_pts_0_quad[] =
{
  { -1.0, -1.0, 0.0 },
  {  1.0, -1.0, 0.0 },
  {  1.0,  1.0, 0.0 },
  { -1.0,  1.0, 0.0 }
};

static double3 lin_pts_1_tri[12] =
{
  {  0.0, -1.0, 0.0 }, // 0
  {  0.0,  0.0, 0.0 }, // 1
  { -1.0,  0.0, 0.0 }, // 2
  { -0.5, -1.0, 0.0 }, // 3
  { -0.5, -0.5, 0.0 }, // 4
  { -1.0, -0.5, 0.0 }, // 5
  {  0.5, -1.0, 0.0 }, // 6
  {  0.5, -0.5, 0.0 }, // 7
  {  0.0, -0.5, 0.0 }, // 8
  { -0.5,  0.0, 0.0 }, // 9
  { -0.5,  0.5, 0.0 }, // 10
  { -1.0,  0.5, 0.0 }  // 11
};

static double3 lin_pts_1_quad[21] =
{
  {  0.0, -1.0, 0.0 }, // 0
  {  1.0,  0.0, 0.0 }, // 1
  {  0.0,  1.0, 0.0 }, // 2
  { -1.0,  0.0, 0.0 }, // 3
  {  0.0,  0.0, 0.0 }, // 4
  { -0.5, -1.0, 0.0 }, // 5
  {  0.0, -0.5, 0.0 }, // 6
  { -0.5,  0.0, 0.0 }, // 7
  { -1.0, -0.5, 0.0 }, // 8
  { -0.5, -0.5, 0.0 }, // 9
  {  0.5, -1.0, 0.0 }, // 10
  {  1.0, -0.5, 0.0 }, // 11
  {  0.5,  0.0, 0.0 }, // 12
  {  0.5, -0.5, 0.0 }, // 13
  {  1.0,  0.5, 0.0 }, // 14
  {  0.5,  1.0, 0.0 }, // 15
  {  0.0,  0.5, 0.0 }, // 16
  {  0.5,  0.5, 0.0 }, // 17
  { -0.5,  1.0, 0.0 }, // 18
  { -1.0,  0.5, 0.0 }, // 19
  { -0.5,  0.5, 0.0 }  // 20
};

int quad_indices[9][5] =
{
  { 0, 1, 2, 3, 4 },
  { 5, 6, 7, 8, 9 }, { 10, 11, 12, 6, 13 },
  { 12, 14, 15, 16, 17 }, { 7, 16, 18, 19, 20 },
  { 0, 11, 4, 8, 6 }, { 4, 14, 2, 19, 16 },
  { 5, 4, 18, 3, 7 }, { 10, 1, 15, 4, 12 }
};

int tri_indices[5][3] =
{
  { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 }, { 9, 10, 11 }, { 9, 4, 8 }
};

int lin_np_tri[2]   = { 3, 12 };
int lin_np_quad[2]  = { 4, 21 };
static int* lin_np[2] = { lin_np_tri, lin_np_quad };

static double3*  lin_tables_tri[2]  = { lin_pts_0_tri, lin_pts_1_tri };
static double3*  lin_tables_quad[2] = { lin_pts_0_quad, lin_pts_1_quad };
static double3** lin_tables[2]      = { lin_tables_tri, lin_tables_quad };


class Quad2DLin : public Quad2D
{
public:

  Quad2DLin()
  {
    mode = HERMES_MODE_TRIANGLE;
    max_order[0]  = max_order[1]  = 1;
    num_tables[0] = num_tables[1] = 2;
    tables = lin_tables;
    np = lin_np;
  };

  virtual void dummy_fn() {}

} quad_lin;



//// vertices and triangles ////////////////////////////////////////////////////////////////////////

Linearizer::Linearizer()
{
  nv = nt = cv = ct = ce = 0;
  verts = NULL;
  tris = NULL;
  edges = NULL;

  pthread_mutexattr_t attr;
  pthread_mutexattr_init(&attr);
  pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
  pthread_mutex_init(&data_mutex, &attr);
  pthread_mutexattr_destroy(&attr);
}


int Linearizer::get_vertex(int p1, int p2, double x, double y, double value)
{
  // search for an existing vertex
  if (p1 > p2) std::swap(p1, p2);
  int index = hash(p1, p2);
  int i = hash_table[index];
  while (i >= 0)
  {
    if (info[i][0] == p1 && info[i][1] == p2 &&
       (value == verts[i][2] || fabs(value - verts[i][2]) < max*1e-4)) return i;
    // note that we won't return a vertex with a different value than the required one;
    // this takes care for discontinuities in the solution, where more vertices
    // with different values will be created
    i = info[i][2];
  }

  // if not found, create a new one
  i = add_vertex();
  verts[i][0] = x;
  verts[i][1] = y;
  verts[i][2] = value;
  info[i][0] = p1;
  info[i][1] = p2;
  info[i][2] = hash_table[index];
  hash_table[index] = i;
  return i;
}


int Linearizer::get_top_vertex(int id, double value)
{
  if (fabs(value - verts[id][2]) < max*1e-24) return id;
  return get_vertex(-rand(), -rand(), verts[id][0], verts[id][1], value);
}


int Linearizer::peek_vertex(int p1, int p2)
{
  // search for a vertex with parents p1, p2
  if (p1 > p2) std::swap(p1, p2);
  int index = hash(p1, p2);
  int i = hash_table[index];
  while (i >= 0)
  {
    if (info[i][0] == p1 && info[i][1] == p2) return i;
    i = info[i][2];
  }
  return -1;
}


void Linearizer::add_triangle(int iv0, int iv1, int iv2)
{
  int index;
  if (del_slot >= 0) // reuse a slot after a deleted triangle
  {
    index = del_slot;
    del_slot = -1;
  }
  else
  {
    if (nt >= ct) {
      tris = (int3*) realloc(tris, sizeof(int3) * (ct = ct * 2));
      verbose("Linearizer::add_triangle(): realloc to %d", ct);
    }
    index = nt++;
  }
  tris[index][0] = iv0;
  tris[index][1] = iv1;
  tris[index][2] = iv2;
}


void Linearizer::print_hash_stats()
{
  const int nh = 10;
  int i, hist[nh] = { 0, };
  for (i = 0; i <= mask; i++)
  {
    int n = 0, j = hash_table[i];
    while (j >= 0 && n < nh-1) { n++; j = info[j][2]; }
    hist[n]++;
  }
  printf("Linearizer: hash histogram: (%d) ", hist[0]);
  for (i = 1; i < nh; i++)
    printf("%d ", hist[i]);
  printf("\n");
}


//// process_triangle & process_quad ///////////////////////////////////////////////////////////////

#ifndef H2D_COMPLEX
  #define getval(i) (val[i])
  #define realpart(x) (x)
#else
  #define getval(i) (val[i].real())
  #define realpart(x) (x.real())
#endif


void Linearizer::process_triangle(int iv0, int iv1, int iv2, int level,
                                  scalar* val, double* phx, double* phy, int* idx)
{
  double midval[3][3];

  if (level < LIN_MAX_LEVEL)
  {
    int i;
    if (!(level & 1))
    {
      // obtain solution values
      sln->set_quad_order(1, item);
      val = sln->get_values(ia, ib);
      if (auto_max)
        for (i = 0; i < lin_np_tri[1]; i++) {
          double v = getval(i);
          if (finite(v) && fabs(v) > max) max = fabs(v);
        }

      // obtain physical element coordinates
      if (curved || disp)
      {
        RefMap* refmap = sln->get_refmap();
        phx = refmap->get_phys_x(1);
        phy = refmap->get_phys_y(1);

        if (disp)
        {
          xdisp->force_transform(sln);
          ydisp->force_transform(sln);
          xdisp->set_quad_order(1, H2D_FN_VAL);
          ydisp->set_quad_order(1, H2D_FN_VAL);
          scalar* dx = xdisp->get_fn_values();
          scalar* dy = ydisp->get_fn_values();
          for (i = 0; i < lin_np_tri[1]; i++) {
            phx[i] += dmult*realpart(dx[i]);
            phy[i] += dmult*realpart(dy[i]);
          }
        }
      }
      idx = tri_indices[0];
    }

    // obtain linearized values and coordinates at the midpoints
    for (i = 0; i < 3; i++)
    {
      midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
      midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
      midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
    };

    // determine whether or not to split the element
    bool split;
    if (eps >= 1.0)
    {
      // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
      split = (level < eps);
    }
    else
    {
      if (!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max && fabs(verts[iv2][2]) > max)
      {
        // do not split if the whole triangle is above the specified maximum value
        split = false;
      }
      else
      {
        // calculate the approximate error of linearizing the normalized solution
        double err = fabs(getval(idx[0]) - midval[2][0]) +
                     fabs(getval(idx[1]) - midval[2][1]) +
                     fabs(getval(idx[2]) - midval[2][2]);
        split = !finite(err) || err > max*3*eps;
      }

      // do the same for the curvature
      if (!split && (curved || disp))
      {
        for (i = 0; i < 3; i++)
          if (sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1.5e-3))
            { split = true; break; }
      }

      // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
      if (level == 0 && !split)
      {
        split = (fabs(getval(8) - 0.5*(midval[2][0] + midval[2][1])) +
                 fabs(getval(9) - 0.5*(midval[2][1] + midval[2][2])) +
                 fabs(getval(4) - 0.5*(midval[2][2] + midval[2][0]))) > max*3*eps;
      }
    }

    // split the triangle if the error is too large, otherwise produce a linear triangle
    if (split)
    {
      if (curved || disp)
        for (i = 0; i < 3; i++) {
          midval[0][i] = phx[idx[i]];
          midval[1][i] = phy[idx[i]];
        }

      // obtain mid-edge vertices
      int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], getval(idx[0]));
      int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], getval(idx[1]));
      int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], getval(idx[2]));

      // recur to sub-elements
      sln->push_transform(0);  process_triangle(iv0, mid0, mid2,  level+1, val, phx, phy, tri_indices[1]);  sln->pop_transform();
      sln->push_transform(1);  process_triangle(mid0, iv1, mid1,  level+1, val, phx, phy, tri_indices[2]);  sln->pop_transform();
      sln->push_transform(2);  process_triangle(mid2, mid1, iv2,  level+1, val, phx, phy, tri_indices[3]);  sln->pop_transform();
      sln->push_transform(3);  process_triangle(mid1, mid2, mid0, level+1, val, phx, phy, tri_indices[4]);  sln->pop_transform();
      return;
    }
  }

  // no splitting: output a linear triangle
  add_triangle(iv0, iv1, iv2);
}


void Linearizer::process_quad(int iv0, int iv1, int iv2, int iv3, int level,
                              scalar* val, double* phx, double* phy, int* idx)
{
  double midval[3][5];

  // try not to split through the vertex with the largest value
  int a = (verts[iv0][2] > verts[iv1][2]) ? iv0 : iv1;
  int b = (verts[iv2][2] > verts[iv3][2]) ? iv2 : iv3;
  a = (verts[a][2] > verts[b][2]) ? a : b;
  int flip = (a == iv1 || a == iv3) ? 1 : 0;

  if (level < LIN_MAX_LEVEL)
  {
    int i;
    if (!(level & 1)) // this is an optimization: do the following only every other time
    {
      // obtain solution values
      sln->set_quad_order(1, item);
      val = sln->get_values(ia, ib);
      if (auto_max)
        for (i = 0; i < lin_np_quad[1]; i++) {
          double v = getval(i);
          if (finite(v) && fabs(v) > max) max = fabs(v);
        }

      // obtain physical element coordinates
      if (curved || disp)
      {
        RefMap* refmap = sln->get_refmap();
        phx = refmap->get_phys_x(1);
        phy = refmap->get_phys_y(1);

        if (disp)
        {
          xdisp->force_transform(sln);
          ydisp->force_transform(sln);
          xdisp->set_quad_order(1, H2D_FN_VAL);
          ydisp->set_quad_order(1, H2D_FN_VAL);
          scalar* dx = xdisp->get_fn_values();
          scalar* dy = ydisp->get_fn_values();
          for (i = 0; i < lin_np_quad[1]; i++) {
            phx[i] += dmult*realpart(dx[i]);
            phy[i] += dmult*realpart(dy[i]);
          }
        }
      }
      idx = quad_indices[0];
    }

    // obtain linearized values and coordinates at the midpoints
    for (i = 0; i < 3; i++)
    {
      midval[i][0] = (verts[iv0][i] + verts[iv1][i]) * 0.5;
      midval[i][1] = (verts[iv1][i] + verts[iv2][i]) * 0.5;
      midval[i][2] = (verts[iv2][i] + verts[iv3][i]) * 0.5;
      midval[i][3] = (verts[iv3][i] + verts[iv0][i]) * 0.5;
      midval[i][4] = (midval[i][0]  + midval[i][2])  * 0.5;
    };

    // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
    midval[2][4] = flip ? (verts[iv0][2] + verts[iv2][2]) * 0.5 : (verts[iv1][2] + verts[iv3][2]) * 0.5;


    // determine whether or not to split the element
    int split;
    if (eps >= 1.0)
    {
      // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
      split = (level < eps) ? 3 : 0;
    }
    else
    {
      if (!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max
                         && fabs(verts[iv2][2]) > max && fabs(verts[iv3][2]) > max)
      {
        // do not split if the whole quad is above the specified maximum value
        split = 0;
      }
      else
      {
        // calculate the approximate error of linearizing the normalized solution
        double herr = fabs(getval(idx[1]) - midval[2][1]) + fabs(getval(idx[3]) - midval[2][3]);
        double verr = fabs(getval(idx[0]) - midval[2][0]) + fabs(getval(idx[2]) - midval[2][2]);
        double err  = fabs(getval(idx[4]) - midval[2][4]) + herr + verr;
        split = (!finite(err) || err > max*4*eps) ? 3 : 0;

        // decide whether to split horizontally or vertically only
        if (level > 0 && split) {
          if (herr > 5*verr)
            split = 1; // h-split
          else if (verr > 5*herr)
            split = 2; // v-split
        }
      }

      // also decide whether to split because of the curvature
      if (split != 3 && (curved || disp))
      {
        double cm2 = sqr(cmax*5e-4);
        if (sqr(phx[idx[1]] - midval[0][1]) + sqr(phy[idx[1]] - midval[1][1]) > cm2 ||
            sqr(phx[idx[3]] - midval[0][3]) + sqr(phy[idx[3]] - midval[1][3]) > cm2) split |= 1;
        if (sqr(phx[idx[0]] - midval[0][0]) + sqr(phy[idx[0]] - midval[1][0]) > cm2 ||
            sqr(phx[idx[2]] - midval[0][2]) + sqr(phy[idx[2]] - midval[1][2]) > cm2) split |= 2;

        /*for (i = 0; i < 5; i++)
          if (sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1e-3))
            { split = 1; break; }*/
      }

      // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
      if (level == 0 && !split)
      {
        split = ((fabs(getval(13) - 0.5*(midval[2][0] + midval[2][1])) +
                  fabs(getval(17) - 0.5*(midval[2][1] + midval[2][2])) +
                  fabs(getval(20) - 0.5*(midval[2][2] + midval[2][3])) +
                  fabs(getval(9)  - 0.5*(midval[2][3] + midval[2][0]))) > max*4*eps) ? 3 : 0;
      }
    }

    // split the quad if the error is too large, otherwise produce two linear triangles
    if (split)
    {
      if (curved || disp)
        for (i = 0; i < 5; i++) {
          midval[0][i] = phx[idx[i]];
          midval[1][i] = phy[idx[i]];
        }

      // obtain mid-edge and mid-element vertices
      int mid0, mid1, mid2, mid3, mid4;
      if (split != 1) mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], getval(idx[0]));
      if (split != 2) mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], getval(idx[1]));
      if (split != 1) mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], getval(idx[2]));
      if (split != 2) mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], getval(idx[3]));
      if (split == 3) mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], getval(idx[4]));

      // recur to sub-elements
      if (split == 3)
      {
        sln->push_transform(0);  process_quad(iv0, mid0, mid4, mid3, level+1, val, phx, phy, quad_indices[1]);  sln->pop_transform();
        sln->push_transform(1);  process_quad(mid0, iv1, mid1, mid4, level+1, val, phx, phy, quad_indices[2]);  sln->pop_transform();
        sln->push_transform(2);  process_quad(mid4, mid1, iv2, mid2, level+1, val, phx, phy, quad_indices[3]);  sln->pop_transform();
        sln->push_transform(3);  process_quad(mid3, mid4, mid2, iv3, level+1, val, phx, phy, quad_indices[4]);  sln->pop_transform();
      }
      else if (split == 1) // h-split
      {
        sln->push_transform(4);  process_quad(iv0, iv1, mid1, mid3, level+1, val, phx, phy, quad_indices[5]);  sln->pop_transform();
        sln->push_transform(5);  process_quad(mid3, mid1, iv2, iv3, level+1, val, phx, phy, quad_indices[6]);  sln->pop_transform();
      }
      else // v-split
      {
        sln->push_transform(6);  process_quad(iv0, mid0, mid2, iv3, level+1, val, phx, phy, quad_indices[7]);  sln->pop_transform();
        sln->push_transform(7);  process_quad(mid0, iv1, iv2, mid2, level+1, val, phx, phy, quad_indices[8]);  sln->pop_transform();
      }
      return;
    }
  }

  // output two linear triangles,
  if (!flip)
  {
    add_triangle(iv3, iv0, iv1);
    add_triangle(iv1, iv2, iv3);
  }
  else
  {
    add_triangle(iv0, iv1, iv2);
    add_triangle(iv2, iv3, iv0);
  }
}


void Linearizer::process_edge(int iv1, int iv2, int marker)
{
  int mid = peek_vertex(iv1, iv2);
  if (mid != -1)
  {
    process_edge(iv1, mid, marker);
    process_edge(mid, iv2, marker);
  }
  else
    add_edge(iv1, iv2, marker);
}


//// regularization ////////////////////////////////////////////////////////////////////////////////

void Linearizer::regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2)
{
  // count the number of hanging mid-edge vertices
  int n = 0;
  if (mid0 >= 0) n++;
  if (mid1 >= 0) n++;
  if (mid2 >= 0) n++;

  if (n == 3)
  {
    // three hanging vertices: split into four triangles
    regularize_triangle(iv0, mid0, mid2, peek_vertex(iv0, mid0), -1, peek_vertex(mid2, iv0));
    regularize_triangle(mid0, iv1, mid1, peek_vertex(mid0, iv1), peek_vertex(iv1, mid1), -1);
    regularize_triangle(mid2, mid1, iv2, -1, peek_vertex(mid1, iv2), peek_vertex(iv2, mid2));
    regularize_triangle(mid0, mid1, mid2, -1, -1, -1);
  }
  else if (n == 2)
  {
    // two hanging vertices: split into three triangles
    if (mid0 < 0)
    {
      regularize_triangle(iv0, iv1, mid1, peek_vertex(iv0, iv1), peek_vertex(iv1, mid1), -1);
      regularize_triangle(mid2, iv0, mid1, peek_vertex(mid2, iv0), -1, -1);
      regularize_triangle(mid2, mid1, iv2, -1, peek_vertex(mid1, iv2), peek_vertex(iv2, mid2));
    }
    else if (mid1 < 0)
    {
      regularize_triangle(iv1, iv2, mid2, peek_vertex(iv1, iv2), peek_vertex(iv2, mid2), -1);
      regularize_triangle(mid0, iv1, mid2, peek_vertex(mid0, iv1), -1, -1);
      regularize_triangle(mid0, mid2, iv0, -1, peek_vertex(mid2, iv0), peek_vertex(iv0, mid0));
    }
    else
    {
      regularize_triangle(iv2, iv0, mid0, peek_vertex(iv2, iv0), peek_vertex(iv0, mid0), -1);
      regularize_triangle(mid1, iv2, mid0, peek_vertex(mid1, iv2), -1, -1);
      regularize_triangle(mid1, mid0, iv1, -1, peek_vertex(mid0, iv1), peek_vertex(iv1, mid1));
    }
  }
  else if (n == 1)
  {
    // one hanging vertex: split into two triangles
    if (mid0 >= 0)
    {
      regularize_triangle(iv0, mid0, iv2, peek_vertex(iv0, mid0), -1, peek_vertex(iv2, iv0));
      regularize_triangle(mid0, iv1, iv2, peek_vertex(mid0, iv1), peek_vertex(iv1, iv2), -1);
    }
    else if (mid1 >= 0)
    {
      regularize_triangle(iv1, mid1, iv0, peek_vertex(iv1, mid1), -1, peek_vertex(iv0, iv1));
      regularize_triangle(mid1, iv2, iv0, peek_vertex(mid1, iv2), peek_vertex(iv2, iv0), -1);
    }
    else
    {
      regularize_triangle(iv2, mid2, iv1, peek_vertex(iv2, mid2), -1, peek_vertex(iv1, iv2));
      regularize_triangle(mid2, iv0, iv1, peek_vertex(mid2, iv0), peek_vertex(iv0, iv1), -1);
    }
  }
  else
  {
    // no hanging vertices: produce a single triangle
    add_triangle(iv0, iv1, iv2);
  }
}


void Linearizer::find_min_max()
{
  // find min & max vertex values
  min_val =  1e100;
  max_val = -1e100;
  for (int i = 0; i < nv; i++)
  {
    if (finite(verts[i][2]) && verts[i][2] < min_val) min_val = verts[i][2];
    if (finite(verts[i][2]) && verts[i][2] > max_val) max_val = verts[i][2];
  }
}


//// process_solution //////////////////////////////////////////////////////////////////////////////

void Linearizer::process_solution(MeshFunction* sln, int item, double eps, double max_abs,
                                  MeshFunction* xdisp, MeshFunction* ydisp, double dmult)
{
  // sanity check
  if (sln == NULL) error("Solution is NULL in Linearizer:process_solution().");

  lock_data();
  TimePeriod time_period;

  // initialization
  this->sln = sln;
  this->item = item;
  this->eps = eps;
  this->xdisp = xdisp;
  this->ydisp = ydisp;
  this->dmult = dmult;
  nv = nt = ne = 0;
  del_slot = -1;

  if (!item) error("Parameter 'item' cannot be zero.");
  get_gv_a_b(item, ia, ib);
  if (ib >= 6) error("Invalid value of parameter 'item'.");

  disp = (xdisp != NULL || ydisp != NULL);
  if (disp && (xdisp == NULL || ydisp == NULL))
    error("Both displacement components must be supplied.");

  // estimate the required number of vertices and triangles
  Mesh* mesh = sln->get_mesh();
  if (mesh == NULL) {
    warn("Have you used Solution::set_coeff_vector() ?");
    error("Mesh is NULL in Linearizer:process_solution().");
  }
  int nn = mesh->get_num_elements();
  int ev = std::max(32 * nn, 10000);  // todo: check this
  int et = std::max(64 * nn, 20000);
  int ee = std::max(24 * nn, 7500);

  // check that displacement meshes are the same
  if (disp)
  {
    unsigned seq1 = mesh->get_seq();
    unsigned seq2 = xdisp->get_mesh()->get_seq();
    unsigned seq3 = ydisp->get_mesh()->get_seq();
    if (seq1 != seq2 || seq1 != seq3)
      error("Displacements must be defined on the same mesh as the solution.");
  }

  // reuse or allocate vertex, triangle and edge arrays
  lin_init_array(verts, double3, cv, ev);
  lin_init_array(tris, int3, ct, et);
  lin_init_array(edges, int3, ce, ee);
  info = (int4*) malloc(sizeof(int4) * cv);

  // initialize the hash table
  int size = 0x2000;
  while (size*2 < cv) size *= 2;
  hash_table = (int*) malloc(sizeof(int) * size);
  memset(hash_table, 0xff, sizeof(int) * size);
  mask = size-1;

  // select the linearization quadrature
  Quad2D *old_quad, *old_quad_x, *old_quad_y;
  old_quad = sln->get_quad_2d();
  sln->set_quad_2d(&quad_lin);
  if (disp) { old_quad_x = xdisp->get_quad_2d();
              old_quad_y = ydisp->get_quad_2d();
              xdisp->set_quad_2d(&quad_lin);
              ydisp->set_quad_2d(&quad_lin); }

  // create all top-level vertices (corresponding to vertex nodes), with
  // all parent-son relations preserved; this is necessary for regularization to
  // work on irregular meshes
  nn = mesh->get_max_node_id();
  int* id2id = new int[nn];
  memset(id2id, 0xff, sizeof(int) * nn);
  bool finished;
  do
  {
    finished = true;
    Node* node;
    for_all_vertex_nodes(node, mesh)
    {
      if (id2id[node->id] < 0 && node->ref != TOP_LEVEL_REF) {
        if (node->p1 < 0)
          id2id[node->id] = get_vertex(node->id, node->id, node->x, node->y, 0);
        else if (id2id[node->p1] >= 0 && id2id[node->p2] >= 0)
          id2id[node->id] = get_vertex(id2id[node->p1], id2id[node->p2], node->x, node->y, 0);
        else
          finished = false;
      }
    }
  }
  while (!finished);

  auto_max = (max_abs < 0.0);
  max = auto_max ? 0.0 : max_abs;

  // obtain the solution in vertices, estimate the maximum solution value
  Element* e;
  for_all_active_elements(e, mesh)
  {
    sln->set_active_element(e);
    sln->set_quad_order(0, item);
    scalar* val = sln->get_values(ia, ib);
    if (val == NULL) error("Item not defined in the solution.");

    scalar *dx, *dy;
    if (disp)
    {
      xdisp->set_active_element(e);
      ydisp->set_active_element(e);
      xdisp->set_quad_order(0, H2D_FN_VAL);
      ydisp->set_quad_order(0, H2D_FN_VAL);
      dx = xdisp->get_fn_values();
      dy = ydisp->get_fn_values();
    }

    for (unsigned int i = 0; i < e->nvert; i++)
    {
      double f = getval(i);
      if (auto_max && finite(f) && fabs(f) > max) max = fabs(f);
      int id = id2id[e->vn[i]->id];
      verts[id][2] = f;

      if (disp)
      {
        verts[id][0] = e->vn[i]->x + dmult*realpart(dx[i]);
        verts[id][1] = e->vn[i]->y + dmult*realpart(dy[i]);
      }
    }
  }

  // process all elements of the mesh
  for_all_active_elements(e, mesh)
  {
    sln->set_active_element(e);
    sln->set_quad_order(0, item);
    scalar* val = sln->get_values(ia, ib);
    if (disp)
    {
      xdisp->set_active_element(e);
      ydisp->set_active_element(e);
    }

    int iv[4];
    for (unsigned int i = 0; i < e->nvert; i++)
      iv[i] = get_top_vertex(id2id[e->vn[i]->id], getval(i));

    // we won't bother calculating physical coordinates from the refmap if this is not a curved element
    curved = e->is_curved();
    cmax = e->get_diameter();

    // recur to sub-elements
    if (e->is_triangle())
      process_triangle(iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL);
    else
      process_quad(iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL);

    for (unsigned int i = 0; i < e->nvert; i++)
      process_edge(iv[i], iv[e->next_vert(i)], e->en[i]->marker);
  }

  delete [] id2id;

  // regularize the linear mesh
  int num = nt;
  for (int i = 0; i < num; i++)
  {
    int iv0 = tris[i][0], iv1 = tris[i][1], iv2 = tris[i][2];
    int mid0 = peek_vertex(iv0, iv1);
    int mid1 = peek_vertex(iv1, iv2);
    int mid2 = peek_vertex(iv2, iv0);
    if (mid0 >= 0 || mid1 >= 0 || mid2 >= 0)
    {
      del_triangle(i);
      regularize_triangle(iv0, iv1, iv2, mid0, mid1, mid2);
    }
  }

  find_min_max();
  //verbose("Linearizer: %d verts, %d tris in %0.3g sec", nv, nt, time_period.tick().last());
  //if (verbose_mode) print_hash_stats();
  unlock_data();

  // select old quadratrues
  sln->set_quad_2d(old_quad);
  if (disp) { xdisp->set_quad_2d(old_quad_x);
              ydisp->set_quad_2d(old_quad_y); }

  // clean up
  ::free(hash_table);
  ::free(info);
}


void Linearizer::free()
{
  lin_free_array(verts, nv, cv);
  lin_free_array(tris, nt, ct);
  lin_free_array(edges, ne, ce);
}


Linearizer::~Linearizer()
{
  free();
  pthread_mutex_destroy(&data_mutex);
}


//// save & load ///////////////////////////////////////////////////////////////////////////////////

void Linearizer::save_data(const char* filename)
{
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open %s for writing.", filename);
  lock_data();

  if (fwrite("H2DL\001\000\000\000", 1, 8, f) != 8 ||
      fwrite(&nv, sizeof(int), 1, f) != 1 ||
      fwrite(verts, sizeof(double3), nv, f) != (unsigned) nv ||
      fwrite(&nt, sizeof(int), 1, f) != 1 ||
      fwrite(tris, sizeof(int3), nt, f) != (unsigned) nt ||
      fwrite(&ne, sizeof(int), 1, f) != 1 ||
      fwrite(edges, sizeof(int3), ne, f) != (unsigned) ne)
  {
    error("Error writing data to %s", filename);
  }

  unlock_data();
  fclose(f);
}


void Linearizer::load_data(const char* filename)
{
  FILE* f = fopen(filename, "rb");
  if (f == NULL) error("Could not open %s for reading.", filename);
  lock_data();

  struct { char magic[4]; int ver; } hdr;
  if (fread(&hdr, sizeof(hdr), 1, f) != 1)
    error("Error reading %s", filename);

  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'L')
    error("File %s is not a Hermes2D Linearizer file.", filename);
  if (hdr.ver > 1)
    error("File %s -- unsupported file version.", filename);

  #define read_array(array, type, n, c, what) \
    if (fread(&n, sizeof(int), 1, f) != 1) \
      error("Error reading the number of " what " from %s", filename); \
    lin_init_array(array, type, c, n); \
    if (fread(array, sizeof(type), n, f) != (unsigned) n) \
      error("Error reading " what " from %s", filename);

  read_array(verts, double3, nv, cv, "vertices");
  read_array(tris,  int3,    nt, ct, "triangles");
  read_array(edges, int3,    ne, ce, "edges");

  find_min_max();
  unlock_data();
  fclose(f);
}

//// others ///////////////////////////////////////////////////////////////////////////////////

void Linearizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const {
  assert_msg(verts != NULL, "Cannot calculate AABB from NULL vertices");
  calc_aabb(&verts[0][0], &verts[0][1], sizeof(double3), nv, min_x, max_x, min_y, max_y);
}

void Linearizer::calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y) {
  *min_x = *max_x = *x;
  *min_y = *max_y = *y;

  uint8_t* ptr_x = (uint8_t*)x;
  uint8_t* ptr_y = (uint8_t*)y;
  for(int i = 0; i < num; i++, ptr_x += stride, ptr_y += stride) {
    *min_x = std::min(*min_x, *((double*)ptr_x));
    *min_y = std::min(*min_y, *((double*)ptr_y));
    *max_x = std::max(*max_x, *((double*)ptr_x));
    *max_y = std::max(*max_y, *((double*)ptr_y));
  }
}
