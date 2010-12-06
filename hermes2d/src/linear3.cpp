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
#include "traverse.h"


extern int tri_indices[5][3];
extern int quad_indices[9][5];

extern int lin_np_tri[2];
extern int lin_np_quad[2];

class Quad2DLin;
extern Quad2DLin quad_lin;


//// vertices and triangles ////////////////////////////////////////////////////////////////////////

Vectorizer::Vectorizer()
          : Linearizer()
{
  verts = NULL;
  dashes = NULL;
  cd = 0;
}


int Vectorizer::create_vertex(double x, double y, double xvalue, double yvalue)
{
  int i = add_vertex();
  verts[i][0] = x;
  verts[i][1] = y;
  verts[i][2] = xvalue;
  verts[i][3] = yvalue;
  return i;
}


int Vectorizer::get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue)
{
  // search for an existing vertex
  if (p1 > p2) std::swap(p1, p2);
  int index = hash(p1, p2);
  int i = hash_table[index];
  while (i >= 0)
  {
    if (info[i][0] == p1 && info[i][1] == p2 &&
       (fabs((xvalue - verts[i][2]) / xvalue) < 1e-4 && fabs((yvalue - verts[i][3]) / yvalue) < 1e-4)) return i;
    // note that we won't return a vertex with a different value than the required one;
    // this takes care for discontinuities in the solution, where more vertices
    // with different values will be created
    i = info[i][2];
  }

  // if not found, create a new one
  i = add_vertex();
  verts[i][0] = x;
  verts[i][1] = y;
  verts[i][2] = xvalue;
  verts[i][3] = yvalue;
  info[i][0] = p1;
  info[i][1] = p2;
  info[i][2] = hash_table[index];
  hash_table[index] = i;
  return i;
}


//// process_triangle & process_quad ///////////////////////////////////////////////////////////////

#ifndef H2D_COMPLEX
  #define getvalx(i) (xval[i])
  #define getvaly(i) (yval[i])
  #define getmag(i) (sqrt(sqr(xval[i]) + sqr(yval[i])))
#else
  #define getvalx(i) (xval[i].real())
  #define getvaly(i) (yval[i].real())
  #define getmag(i) (sqrt(sqr(xval[i].real()) + sqr(yval[i].real())))
#endif
#define midmag(i) (sqrt(sqr(midval[2][i]) + sqr(midval[3][i])))
#define magvert(i) (sqr(verts[i][2]) + sqr(verts[i][3]))


void Vectorizer::process_triangle(int iv0, int iv1, int iv2, int level,
                                  scalar* xval, scalar* yval, double* phx, double* phy, int* idx)
{
  if (level < LIN_MAX_LEVEL)
  {
    int i;
    if (!(level & 1))
    {
      // obtain solution values and physical element coordinates
      xsln->set_quad_order(1, xitem);
      ysln->set_quad_order(1, yitem);
      xval = xsln->get_values(xia, xib);
      yval = ysln->get_values(yia, yib);
      for (i = 0; i < lin_np_tri[1]; i++) {
        double m = getmag(i);
        if (finite(m) && fabs(m) > max) max = fabs(m);
      }

      if (curved)
      {
        RefMap* refmap = xsln->get_refmap();
        phx = refmap->get_phys_x(1);
        phy = refmap->get_phys_y(1);
      }
      idx = tri_indices[0];
    }

    // obtain linearized values and coordinates at the midpoints
    double midval[4][3];
    for (i = 0; i < 4; i++)
    {
      midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
      midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
      midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
    };

    // determine whether or not to split the element
    bool split;
    if (eps >= 1.0)
    {
      //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
      split = (level < eps);
    }
    else
    {
      // calculate the approximate error of linearizing the normalized solution
      double err = fabs(getmag(idx[0]) - midmag(0)) +
                   fabs(getmag(idx[1]) - midmag(1)) +
                   fabs(getmag(idx[2]) - midmag(2));
      split = !finite(err) || err > max*3*eps;

      // do the same for the curvature
      if (curved && !split)
      {
        double cerr = 0.0, cden = 0.0; // fixme
        for (i = 0; i < 3; i++)
        {
          cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
          cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
        }
        split = (cerr > cden*2.5e-4);
      }

      // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
      if (level == 0 && !split)
      {
        split = (fabs(getmag(8) - 0.5*(midmag(0) + midmag(1))) +
                 fabs(getmag(9) - 0.5*(midmag(1) + midmag(2))) +
                 fabs(getmag(4) - 0.5*(midmag(2) + midmag(0)))) > max*3*eps;
      }
    }

    // split the triangle if the error is too large, otherwise produce a linear triangle
    if (split)
    {
      if (curved)
        for (i = 0; i < 3; i++) {
          midval[0][i] = phx[idx[i]];
          midval[1][i] = phy[idx[i]];
        }

      // obtain mid-edge vertices
      int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], getvalx(idx[0]), getvaly(idx[0]));
      int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], getvalx(idx[1]), getvaly(idx[1]));
      int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], getvalx(idx[2]), getvaly(idx[2]));

      // recur to sub-elements
      push_transform(0);  process_triangle(iv0, mid0, mid2,  level+1, xval, yval, phx, phy, tri_indices[1]);  pop_transform();
      push_transform(1);  process_triangle(mid0, iv1, mid1,  level+1, xval, yval, phx, phy, tri_indices[2]);  pop_transform();
      push_transform(2);  process_triangle(mid2, mid1, iv2,  level+1, xval, yval, phx, phy, tri_indices[3]);  pop_transform();
      push_transform(3);  process_triangle(mid1, mid2, mid0, level+1, xval, yval, phx, phy, tri_indices[4]);  pop_transform();
      return;
    }
  }

  // no splitting: output a linear triangle
  add_triangle(iv0, iv1, iv2);
}


void Vectorizer::process_quad(int iv0, int iv1, int iv2, int iv3, int level,
                              scalar* xval, scalar* yval, double* phx, double* phy, int* idx)
{
  // try not to split through the vertex with the largest value
  int a = (magvert(iv0) > magvert(iv1)) ? iv0 : iv1;
  int b = (magvert(iv2) > magvert(iv3)) ? iv2 : iv3;
  a = (magvert(a) > magvert(b)) ? a : b;

  int flip = (a == iv1 || a == iv3) ? 1 : 0;

  if (level < LIN_MAX_LEVEL)
  {
    int i;
    if (!(level & 1))
    {
      // obtain solution values and physical element coordinates
      xsln->set_quad_order(1, xitem);
      ysln->set_quad_order(1, yitem);
      xval = xsln->get_values(xia, xib);
      yval = ysln->get_values(yia, yib);
      for (i = 0; i < lin_np_quad[1]; i++) {
        double m = getmag(i);
        if (finite(m) && fabs(m) > max) max = fabs(m);
      }

      if (curved)
      {
        RefMap* refmap = xsln->get_refmap();
        phx = refmap->get_phys_x(1);
        phy = refmap->get_phys_y(1);
      }
      idx = quad_indices[0];
    }

    // obtain linearized values and coordinates at the midpoints
    double midval[4][5];
    for (i = 0; i < 4; i++)
    {
      midval[i][0] = (verts[iv0][i] + verts[iv1][i]) * 0.5;
      midval[i][1] = (verts[iv1][i] + verts[iv2][i]) * 0.5;
      midval[i][2] = (verts[iv2][i] + verts[iv3][i]) * 0.5;
      midval[i][3] = (verts[iv3][i] + verts[iv0][i]) * 0.5;
      midval[i][4] = (midval[i][0]  + midval[i][2])  * 0.5;
    };

    // determine whether or not to split the element
    bool split;
    if (eps >= 1.0)
    {
      //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
      split = (level < eps);
    }

    else
    {
      // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
      midval[2][4] = flip ? (verts[iv0][2] + verts[iv2][2]) * 0.5 : (verts[iv1][2] + verts[iv3][2]) * 0.5;
      midval[3][4] = flip ? (verts[iv0][3] + verts[iv2][3]) * 0.5 : (verts[iv1][3] + verts[iv3][3]) * 0.5;

      // calculate the approximate error of linearizing the normalized solution
      double err = fabs(getmag(idx[0]) - midmag(0)) +
                   fabs(getmag(idx[1]) - midmag(1)) +
                   fabs(getmag(idx[2]) - midmag(2)) +
                   fabs(getmag(idx[3]) - midmag(3)) +
                   fabs(getmag(idx[4]) - midmag(4));
      split = !finite(err) || err > max*4*eps;

      // do the same for the curvature
      if (curved && !split)
      {
        double cerr = 0.0, cden = 0.0; // fixme
        for (i = 0; i < 5; i++)
        {
          cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
          cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
        }
        split = (cerr > cden*2.5e-4);
      }

      // do extra tests at level 0, so as not to miss functions with zero error at edge midpoints
      if (level == 0 && !split)
      {
        err = fabs(getmag(13) - 0.5*(midmag(0) + midmag(1))) +
              fabs(getmag(17) - 0.5*(midmag(1) + midmag(2))) +
              fabs(getmag(20) - 0.5*(midmag(2) + midmag(3))) +
              fabs(getmag(9)  - 0.5*(midmag(3) + midmag(0)));
        split = !finite(err) || (err) > max*2*eps; //?
      }
    }

    // split the quad if the error is too large, otherwise produce two linear triangles
    if (split)
    {
      if (curved)
        for (i = 0; i < 5; i++) {
          midval[0][i] = phx[idx[i]];
          midval[1][i] = phy[idx[i]];
        }

      // obtain mid-edge and mid-element vertices
      int mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], getvalx(idx[0]), getvaly(idx[0]));
      int mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], getvalx(idx[1]), getvaly(idx[1]));
      int mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], getvalx(idx[2]), getvaly(idx[2]));
      int mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], getvalx(idx[3]), getvaly(idx[3]));
      int mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], getvalx(idx[4]), getvaly(idx[4]));

      // recur to sub-elements
      push_transform(0);  process_quad(iv0, mid0, mid4, mid3, level+1, xval, yval, phx, phy, quad_indices[1]);  pop_transform();
      push_transform(1);  process_quad(mid0, iv1, mid1, mid4, level+1, xval, yval, phx, phy, quad_indices[2]);  pop_transform();
      push_transform(2);  process_quad(mid4, mid1, iv2, mid2, level+1, xval, yval, phx, phy, quad_indices[3]);  pop_transform();
      push_transform(3);  process_quad(mid3, mid4, mid2, iv3, level+1, xval, yval, phx, phy, quad_indices[4]);  pop_transform();
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


void Vectorizer::process_dash(int iv1, int iv2)
{
  int mid = peek_vertex(iv1, iv2);
  if (mid != -1)
  {
    process_dash(iv1, mid);
    process_dash(mid, iv2);
  }
  else
    add_dash(iv1, iv2);
}


void Vectorizer::find_min_max()
{
  // find min & max vertex values
  min_val =  1e100;
  max_val = -1e100;
  for (int i = 0; i < nv; i++)
  {
    double mag = (verts[i][2]*verts[i][2] + verts[i][3]*verts[i][3]);
    if (finite(mag) && mag < min_val) min_val = mag;
    if (finite(mag) && mag > max_val) max_val = mag;
  }
  max_val = sqrt(max_val);
  min_val = sqrt(min_val);

}


//// process_solution //////////////////////////////////////////////////////////////////////////////

void Vectorizer::process_solution(MeshFunction* xsln, int xitem, MeshFunction* ysln, int yitem, double eps)
{
  // sanity check
  if (xsln == NULL || ysln == NULL) error("One of the solutions is NULL in Vectorizer:process_solution().");


  lock_data();
  TimePeriod cpu_time;

  // initialization
  this->xsln = xsln;
  this->ysln = ysln;
  this->xitem = xitem;
  this->yitem = yitem;
  this->eps = eps;
  nv = nt = ne = nd = 0;
  del_slot = -1;

  Mesh* meshes[2] = { xsln->get_mesh(), ysln->get_mesh() };
  if (meshes[0] == NULL || meshes[1] == NULL) {
    error("One of the meshes is NULL in Vectorizer:process_solution().");
  }

  Transformable* fns[2] = { xsln, ysln };
  Traverse trav;

  // estimate the required number of vertices and triangles
  // (based on the assumption that the linear mesh will be
  // about four-times finer than the original mesh).
  int nn = meshes[0]->get_num_elements() + meshes[1]->get_num_elements();
  int ev = std::max(32 * nn, 10000);
  int et = std::max(64 * nn, 20000);
  int ee = std::max(24 * nn, 7500);
  int ed = ee;

  lin_init_array(verts, double4, cv, ev);
  lin_init_array(tris, int3, ct, et);
  lin_init_array(edges, int3, ce, ee);
  lin_init_array(dashes, int2, cd, ed);

  info = (int4*) malloc(sizeof(int4) * cv);

  // initialize the hash table
  int size = 0x1000;
  while (size*2 < cv) size *= 2;
  hash_table = (int*) malloc(sizeof(int) * size);
  memset(hash_table, 0xff, sizeof(int) * size);
  mask = size-1;


  // select the linearization quadrature
  Quad2D *old_quad_x, *old_quad_y;
  old_quad_x = xsln->get_quad_2d();
  old_quad_y = ysln->get_quad_2d();

  xsln->set_quad_2d((Quad2D*) &quad_lin);
  ysln->set_quad_2d((Quad2D*) &quad_lin);

  if (!xitem) error("Parameter 'xitem' cannot be zero.");
  if (!yitem) error("Parameter 'yitem' cannot be zero.");
  get_gv_a_b(xitem, xia, xib);
  get_gv_a_b(yitem, yia, yib);
  if (xib >= 6) error("Invalid value of paremeter 'xitem'.");
  if (yib >= 6) error("Invalid value of paremeter 'yitem'.");

  max = 1e-10;
  trav.begin(2, meshes, fns);
  Element** e;
  while ((e = trav.get_next_state(NULL, NULL)) != NULL)
  {
    xsln->set_quad_order(0, xitem);
    ysln->set_quad_order(0, yitem);
    scalar* xval = xsln->get_values(xia, xib);
    scalar* yval = ysln->get_values(yia, yib);

    for (unsigned int i = 0; i < e[0]->nvert; i++)
    {
      double fx = getvalx(i);
      double fy = getvaly(i);
      if (fabs(sqrt(fx*fx + fy*fy)) > max) max = fabs(sqrt(fx*fx + fy*fy));
    }
  }
  trav.finish();

  trav.begin(2, meshes, fns);
  // process all elements of the mesh
  while ((e = trav.get_next_state(NULL, NULL)) != NULL)
  {
    xsln->set_quad_order(0, xitem);
    ysln->set_quad_order(0, yitem);
    scalar* xval = xsln->get_values(xia, xib);
    scalar* yval = ysln->get_values(yia, yib);

    double* x = xsln->get_refmap()->get_phys_x(0);
    double* y = ysln->get_refmap()->get_phys_y(0);

    int iv[4];
    for (unsigned int i = 0; i < e[0]->nvert; i++)
    {
      double fx = getvalx(i);
      double fy = getvaly(i);
      iv[i] = create_vertex(x[i], y[i], fx, fy);
    }

    // we won't bother calculating physical coordinates from the refmap if this is not a curved element
    curved = (e[0]->cm != NULL);

    // recur to sub-elements
    if (e[0]->is_triangle())
      process_triangle(iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL, NULL);
    else
      process_quad(iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL, NULL);

    // process edges and dashes (bold line for edge in both meshes, dashed line for edge in one of the meshes)
    Trf* xctm = xsln->get_ctm();
    Trf* yctm = ysln->get_ctm();
    double r[4] = { -1.0, 1.0, 1.0, -1.0 };
    double ref[4][2] = { {-1.0,-1.0}, {1.0,-1.0}, {1.0,1.0}, {-1.0,1.0} };
    for (unsigned int i = 0; i < e[0]->nvert; i++)
    {
      bool bold = false;
      double px = ref[i][0];
      double py = ref[i][1];
      // for odd edges (1, 3) we check x coordinate after ctm transformation, if it's the same (1 or -1) in both meshes => bold
      if (i & 1) {
        if ((xctm->m[0]*px + xctm->t[0] == r[i]) && (yctm->m[0]*px + yctm->t[0] == r[i]))
          bold = true;
      }
      // for even edges (0, 4) we check y coordinate after ctm transformation, if it's the same (-1 or 1) in both meshes => bold
      else {
        if ((xctm->m[1]*py + xctm->t[1] == r[i]) && (yctm->m[1]*py + yctm->t[1] == r[i]))
          bold = true;
      }
      int j = e[0]->next_vert(i);
      // we draw a line only if both edges lies on the boundary or if the line is from left top to right bottom
      if (((e[0]->en[i]->bnd) && (e[1]->en[i]->bnd)) ||
         (verts[iv[i]][1] < verts[iv[j]][1]) ||
         (verts[iv[i]][1] == verts[iv[j]][1] && verts[iv[i]][0] < verts[iv[j]][0]))
      {
        if (bold)
          process_edge(iv[i], iv[j], e[0]->en[i]->marker);
        else
          process_dash(iv[i], iv[j]);
      }
    }
  }
  trav.finish();

  find_min_max();

  verbose("Vectorizer created %d verts and %d tris in %0.3g s", nv, nt, cpu_time.tick().last());
  //if (verbose_mode) print_hash_stats();
  unlock_data();

   // select old quadratrues
  xsln->set_quad_2d(old_quad_x);
  ysln->set_quad_2d(old_quad_y);

  // clean up
  ::free(hash_table);
  ::free(info);

}


//// save & load ///////////////////////////////////////////////////////////////////////////////////

void Vectorizer::save_data(const char* filename)
{
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open %s for writing.", filename);
  lock_data();

  if (fwrite("H2DV\001\000\000\000", 1, 8, f) != 8 ||
      fwrite(&nv, sizeof(int), 1, f) != 1 ||
      fwrite(verts, sizeof(double4), nv, f) != (unsigned) nv ||
      fwrite(&nt, sizeof(int), 1, f) != 1 ||
      fwrite(tris, sizeof(int3), nt, f) != (unsigned) nt ||
      fwrite(&ne, sizeof(int), 1, f) != 1 ||
      fwrite(edges, sizeof(int3), ne, f) != (unsigned) ne ||
      fwrite(&nd, sizeof(int), 1, f) != 1 ||
      fwrite(dashes, sizeof(int2), nd, f) != (unsigned) nd)
  {
    error("Error writing data to %s", filename);
  }

  unlock_data();
  fclose(f);
}


void Vectorizer::load_data(const char* filename)
{
  FILE* f = fopen(filename, "rb");
  if (f == NULL) error("Could not open %s for reading.", filename);
  lock_data();

  struct { char magic[4]; int ver; } hdr;
  if (fread(&hdr, sizeof(hdr), 1, f) != 1)
    error("Error reading %s", filename);

  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'V')
    error("File %s is not a Hermes2D Vectorizer file.", filename);
  if (hdr.ver > 1)
    error("File %s -- unsupported file version.", filename);

  #define read_array(array, type, n, c, what) \
    if (fread(&n, sizeof(int), 1, f) != 1) \
      error("Error reading the number of " what " from %s", filename); \
    lin_init_array(array, type, c, n); \
    if (fread(array, sizeof(type), n, f) != (unsigned) n) \
      error("Error reading " what " from %s", filename);

  read_array(verts, double4, nv, cv, "vertices");
  read_array(tris,  int3,    nt, ct, "triangles");
  read_array(edges, int3,    ne, ce, "edges");
  read_array(dashes, int2,   nd, cd, "dashes");

  find_min_max();
  unlock_data();
  fclose(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////

Vectorizer::~Vectorizer()
{
  lin_free_array(verts, nv, cv);
  // fixme: free dashes?
}

//// others ///////////////////////////////////////////////////////////////////////////////////

void Vectorizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const {
  assert_msg(verts != NULL, "Cannot calculate AABB from NULL vertices");
  calc_aabb(&verts[0][0], &verts[0][1], sizeof(double4), nv, min_x, max_x, min_y, max_y);
}
