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

#include "../h2d_common.h"
#include "shapeset.h"
#include "../../../hermes_common/matrix.h"


/*    numbering of edge intervals: (the variable 'part')

        -+-        -+-         -+-
         |          |        13 |
         |        5 |          -+-
         |          |        12 |
       1 |         -+-         -+-                  finer interval:
         |          |        11 |
         |        4 |          -+-                  p = (p + 1) * 2   (+1)
         |          |        10 |
        -+-        -+-         -+-   ... etc.
         |          |         9 |
         |        3 |          -+-
         |          |         8 |
       0 |         -+-         -+-
         |          |         7 |
         |        2 |          -+-
         |          |         6 |
        -+-        -+-         -+-            */


/// Constrained edge functions are constructed by subtracting the linear part (ie., two
/// vertex functions) from the constraining edge function and expressing the rest as a
/// linear combination of standard edge functions. This function determines the coefficients
/// of such linear combination by forming and solving a simple linear system.
///
double* Shapeset::calculate_constrained_edge_combination(int order, int part, int ori)
{
  /*
  "ebias" is the order of the first edge function, this has to be subtracted
  from the order to get a reasonable numbering of the edge functions, starting
  from 0. For H1 ebias is 2 and for Hcurl it is 0.
  */
  assert((order - ebias) >= 0);
  assert(part >= 0);

  int i, j, n;

  // determine the interval of the edge
  for (n = 2; n <= part; n <<= 1)
    part -= n;

  double n2 = 2.0 / n;
  double hi = -((double) part * n2 - 1.0);
  double lo = -((double) (part + 1) * n2 - 1.0);

  int idx[16];
  ori = ori ? 0 : 1;
  for (i = 0; i <= order; i++)
    idx[i] = get_edge_index(0, ori, i);

  // function values at the endpoints (for subtracting of the linear part in H1)
  double c = 1.0;
  double f_lo = 0.0, f_hi = 0.0;
  if (ebias == 2)
  {
    f_lo = get_value(0, idx[order], lo, -1.0, 0);
    f_hi = get_value(0, idx[order], hi, -1.0, 0);
  }
  else
  {
    // this is for H(curl), where we need to multiply the constrained function
    // so that its vectors are not shortened and fit the large constraining fn.
    c = (hi - lo) / 2.0;
  }

  // fill the matrix of the linear system
  n = order + 1 - ebias;
  int space_type = this->get_space_type();
  int component = (space_type == HERMES_HDIV_SPACE)? 1 : 0;
  double** a = new_matrix<double>(n, n);
  double* b = new double[n];
  for (i = 0; i < n; i++)
  {
    // chebyshev point
    int o = (ebias == 0) ? order + 1 : order;
    double p = cos((i+1) * M_PI / o);
    double r = (p + 1.0) * 0.5;
    double s = 1.0 - r;

    // matrix row
    for (j = 0; j < n; j++)
      a[i][j] = get_value(0, idx[j+ebias], p, -1.0, component);

    // rhs
    b[i] = c * get_value(0, idx[order], lo*s + hi*r, -1.0, component) - f_lo*s - f_hi*r;
  }

  // solve the system
  double d;
  int* iperm = new int[n];
  ludcmp(a, n, iperm, &d);
  lubksb(a, n, iperm, b);

  // cleanup
  delete [] iperm;
  delete [] a;

  return b;
}


/// Returns the coefficients for the linear combination forming a constrained edge function.
/// This function performs the storage (caching) of these coefficients, so that they can be
/// calculated only once.
///
double* Shapeset::get_constrained_edge_combination(int order, int part, int ori, int& nitems)
{
  int index = 2*((max_order + 1 - ebias)*part + (order - ebias)) + ori;

  // allocate/reallocate the array if necessary
  if (comb_table == NULL)
  {
    table_size = 1024;
    while (table_size <= index) table_size *= 2;
    comb_table = (double**) malloc(table_size * sizeof(double*));
    memset(comb_table, 0, table_size * sizeof(double*));
  }
  else if (index >= table_size)
  {
    // adjust table_size to accommodate the required depth
    int old_size = table_size;
    while (index >= table_size) table_size *= 2;

    // reallocate the table
    verbose("Shapeset::get_constrained_edge_combination(): realloc to table_size=%d", table_size);
    comb_table = (double**) realloc(comb_table, table_size * sizeof(double*));
    memset(comb_table + old_size, 0,(table_size - old_size) * sizeof(double*));
  }

  // do we have the required linear combination yet?
  if (comb_table[index] == NULL)
  {
    // no, calculate it
    comb_table[index] = calculate_constrained_edge_combination(order, part, ori);
  }

  nitems = order + 1 - ebias;
  return comb_table[index];
}


/// Releases all cached coefficients.
void Shapeset::free_constrained_edge_combinations()
{
  if (comb_table != NULL)
  {
    for (int i = 0; i < table_size; i++)
      if (comb_table[i] != NULL)
        delete [] comb_table[i];

    free(comb_table);
    comb_table = NULL;
  }
}


#define parse_index \
    int part = (unsigned) index >> 7, \
        order = (index >> 3) & 15, \
        edge = (index >> 1) & 3, \
        ori = index & 1;


/// Constructs the linear combination of edge functions, forming a constrained edge function.
///
double Shapeset::get_constrained_value(int n, int index, double x, double y, int component)
{
  index = -1 - index;
  parse_index;

  int i, nc;
  double sum, *comb = get_constrained_edge_combination(order, part, ori, nc);

  sum = 0.0;
  shape_fn_t* table = shape_table[n][mode][component];
  for (i = 0; i < nc; i++)
    sum += comb[i] * table[get_edge_index(edge, ori, i+ebias)](x, y);

  return sum;
}
