// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include "config.h"

#include <string.h>
#include "mesh_reader_exodusii.h"
#include "mesh.h"
#include <map>

#ifdef WITH_EXODUSII
#include <exodusII.h>
#endif
namespace Hermes
{
  namespace Hermes2D
  {
    MeshReaderExodusII::MeshReaderExodusII()
    {
#ifdef WITH_EXODUSII
#else
      throw Hermes::Exceptions::Exception("hermes2d was not compiled with ExodusII support");
#endif
    }

    MeshReaderExodusII::~MeshReaderExodusII()
    {
    }

    struct Vertex
    {
      double x, y;
    };

    // needed for comparison of vertices inside std::map
    struct VCompare
    {
      bool operator()(Vertex a, Vertex b) const
      {
        if(a.x < b.x) return true;
        else if(a.x > b.x) return false;
        else
        {
          if(a.y < b.y) return true;
          return false;
        }
      }
    };

    bool MeshReaderExodusII::load(const char *file_name, Mesh *mesh)
    {
#ifdef WITH_EXODUSII
      int err;
      int cpu_ws = sizeof(double);    // use float or double
      int io_ws = 8;            // store variables as doubles
      float version;
      int exoid = ex_open(file_name, EX_READ, &cpu_ws, &io_ws, &version);

      // read initialization parameters
      int n_dims, n_nodes, n_elems, n_eblocks, n_nodesets, n_sidesets;
      char title[MAX_LINE_LENGTH + 1];
      err = ex_get_init(exoid, title, &n_dims, &n_nodes, &n_elems, &n_eblocks, &n_nodesets, &n_sidesets);
      if(n_dims != 2)
      {
        throw Hermes::Exceptions::Exception("File '%s' does not contain 2D mesh", file_name);
        return false;
      }

      // load coordinates
      double *x = new double[n_nodes];
      double *y = new double[n_nodes];
      err = ex_get_coord(exoid, x, y, NULL);

      // remove duplicate vertices and build renumbering map
      std::map<Vertex, int, VCompare> vtx_list;        // map for eliminating duplicities
      std::map<int, int> vmap;                // reindexing map
      Hermes::vector<Vertex> vtx_arr;              // vertex array
      int vid = 0;
      for (int i = 0; i < n_nodes; i++)
      {
        int k;
        Vertex v = { x[i], y[i] };
        if(vtx_list.count(v) == 0)
        {
          vtx_arr.push_back(v);
          k = vid++;
          vtx_list[v] = k;
        }
        else
          k = vtx_list[v];

        vmap[i + 1] = k;
      }
      delete [] x;
      delete [] y;

      int n_vtx = vtx_arr.size();
      double2 *vtx = new double2[n_vtx];
      for (int i = 0; i < n_vtx; i++)
      {
        vtx[i][0] = vtx_arr[i].x;
        vtx[i][1] = vtx_arr[i].y;
      }

      int n_tri = 0;    // number of triangles
      int n_quad = 0;    // number of quads

      // get info about element blocks
      int *eid_blocks = new int[n_eblocks];
      err = ex_get_elem_blk_ids(exoid, eid_blocks);
      // go over all element blocks
      for (int i = 0; i < n_eblocks; i++)
      {
        int id = eid_blocks[i];

        // get block info
        char elem_type[MAX_STR_LENGTH + 1];
        int n_elems_in_blk, n_elem_nodes, n_attrs;
        err = ex_get_elem_block(exoid, id, elem_type, &n_elems_in_blk, &n_elem_nodes, &n_attrs);

        if(n_elem_nodes == 3) n_tri += n_elems_in_blk;
        else if(n_elem_nodes == 4) n_quad += n_elems_in_blk;
        else
        {
          delete [] vtx;
          throw Hermes::Exceptions::Exception("Unknown type of element");
          return false;
        }
      }
      int3 *tri = n_tri > 0 ? new int3[n_tri] : NULL;    // triangles
      std::string *tri_markers = n_tri > 0 ? new std::string[n_tri] : NULL;
      int4 *quad = n_quad > 0 ? new int4[n_quad] : NULL;    // quads
      std::string *quad_markers = n_quad > 0 ? new std::string[n_quad] : NULL;

      int n_els = n_tri + n_quad;                // total number of elements
      int **els = n_els > 0 ? new int *[n_els] : NULL;    // elements
      int *el_nv = n_els > 0 ? new int[n_els] : NULL;    // number of vertices for each element

      int it = 0, iq = 0, iel = 0;
      for (int i = 0; i < n_eblocks; i++)
      {
        int id = eid_blocks[i];

        // get block info
        char elem_type[MAX_STR_LENGTH + 1];
        int n_elems_in_blk, n_elem_nodes, n_attrs;
        err = ex_get_elem_block(exoid, id, elem_type, &n_elems_in_blk, &n_elem_nodes, &n_attrs);

        // read connectivity array
        int *connect = new int[n_elem_nodes * n_elems_in_blk];
        err = ex_get_elem_conn(exoid, id, connect);

        // Update the mesh' internal array element_markers_conversion.
        std::ostringstream string_stream;
        string_stream << id;
        std::string el_marker = string_stream.str();

        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->element_markers_conversion.insert_marker(mesh->element_markers_conversion.min_marker_unused, el_marker);

        int ic = 0;
        for (int j = 0; j < n_elems_in_blk; j++)
        {
          el_nv[iel] = n_elem_nodes;
          if(n_elem_nodes == 3)
          {
            tri[it][0] = vmap[connect[ic++]];
            tri[it][1] = vmap[connect[ic++]];
            tri[it][2] = vmap[connect[ic++]];
            tri_markers[it] = el_marker;
            els[iel] = tri[it];
            it++;
          }
          else if(n_elem_nodes == 4)
          {
            quad[iq][0] = vmap[connect[ic++]];
            quad[iq][1] = vmap[connect[ic++]];
            quad[iq][2] = vmap[connect[ic++]];
            quad[iq][3] = vmap[connect[ic++]];
            quad_markers[iq] = el_marker;
            els[iel] = quad[iq];
            iq++;
          }
          else
          {
            delete [] vtx;
            throw Hermes::Exceptions::Exception("Unknown type of element");
            return false;
          }
          iel++;
        }
        delete [] connect;
      }
      delete [] eid_blocks;

      // query number of side sets
      int *sid_blocks = new int[n_sidesets];
      err = ex_get_side_set_ids(exoid, sid_blocks);

      // go over the sidesets
      int n_mark = 0;    // number of markers
      for (int i = 0; i < n_sidesets; i++)
      {
        int sid = sid_blocks[i];
        int n_sides_in_set, n_df_in_set;
        err = ex_get_side_set_param(exoid, sid, &n_sides_in_set, &n_df_in_set);
        n_mark += n_sides_in_set;
      }
      int2 *marks = new int2[n_mark];
      std::string *bnd_markers = new std::string[n_mark];

      int im = 0;
      for (int i = 0; i < n_sidesets; i++)
      {
        int sid = sid_blocks[i];
        int n_sides_in_set, n_df_in_set;
        err = ex_get_side_set_param(exoid, sid, &n_sides_in_set, &n_df_in_set);
        int num_elem_in_set = n_sides_in_set;

        int *elem_list = new int[num_elem_in_set];
        int *side_list = new int[n_sides_in_set];
        err = ex_get_side_set(exoid, sid, elem_list, side_list);

        // Update the mesh' internal array boundary_markers_conversion.
        std::ostringstream string_stream;
        string_stream << sid;
        std::string bnd_marker = string_stream.str();

        // This functions check if the user-supplied marker on this element has been
        // already used, and if not, inserts it in the appropriate structure.
        mesh->boundary_markers_conversion.insert_marker(mesh->boundary_markers_conversion.min_marker_unused, bnd_marker);

        for (int j = 0; j < num_elem_in_set; j++)
        {
          int nv = el_nv[side_list[j] - 1];      // # of vertices of the element
          int vt = side_list[j] - 1;
          marks[im][0] = els[elem_list[j] - 1][vt];
          marks[im][1] = els[elem_list[j] - 1][(vt + 1) % nv];
          bnd_markers[im] = bnd_marker;
          im++;
        }

        delete [] elem_list;
        delete [] side_list;
      }
      delete [] sid_blocks;

      // we are done
      err = ex_close(exoid);

      mesh->create(n_vtx, vtx, n_tri, tri, tri_markers, n_quad, quad, quad_markers, n_mark, marks, bnd_markers);

      // clean-up
      delete [] marks;
      delete [] tri;
      delete [] quad;
      delete [] tri_markers;
      delete [] quad_markers;
      delete [] bnd_markers;
      delete [] vtx;
      delete [] el_nv;
      delete [] els;

      return true;
#else
      return false;
#endif
    }
  }
}