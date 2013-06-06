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

#ifndef __H2D_MESH_UTIL_H
#define __H2D_MESH_UTIL_H

namespace Hermes
{
  namespace Hermes2D
  {
    /*  node and son numbering on a triangle:

    -Triangle to triangles refinement

    vn[2]                                       vn[2]

    *                                           *

    / \                                         / \
    /   \                                       /   \
    /     \                                     /     \
    /       \                                   / son[2]\
    /         \                                 /_________\
    en[2]   /           \   en[1]                 vn[0] *           * vn[1]
    *             *                       vn[1]  *-----------*  vn[0]
    /               \                     vn[2] *  \         /  * vn[2]
    /                 \                         / \  \ son[3]/  / \
    /                   \                       /   \  \     /  /   \
    /                     \                     /     \  \   /  /     \
    /                       \                   / son[0]\  \ /  /son[1] \
    /                         \                 /         \  *  /         \
    *-------------*-------------*               *-----------*   *-----------*
    vn[0]      vn[1] vn[2] vn[0]      vn[1]
    vn[0]           en[0]           vn[1]

    -Triangle to quads refinement

    vn[2]                                     vn[2]

    *                                        *
    / \                                      / \
    /   \                                    /   \
    /     \                                  /     \
    /       \                          vn[3] * son[2]* vn[1]
    /         \                       vn[3] *  \     /  * vn[2]
    en[2]   *           *   en[1]                   / \  \   /  / \
    /             \                         /   \ vn[0] /   \
    /               \                       /     \  *  /     \
    /                 \                     /       \   /       \
    /         *         \                   /   vn[2] * * vn[3]   \
    /                     \                 /          | |          \
    /                       \               /  son[0]   | |  son[1]   \
    /                         \             /            | |            \
    *-------------*-------------*           *-------------* *-------------*
    vn[0]      vn[1]   vn[0]        vn[1]
    vn[0]           en[0]           vn[1]

    node and son numbering on a quad:          refinement '0':

    vn[3]           en[2]           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

    *-------------*-------------*               *------------* *------------*
    |                           |               |            | |            |
    |                           |               |            | |            |
    |                           |               |   son[3]   | |   son[2]   |
    |                           |               |            | |            |
    |                           |               |       vn[1]| |vn[0]       |
    |                           |         vn[0] *------------* *------------* vn[1]
    en[3]  *                           *  en[1]  vn[3] *------------* *------------* vn[2]
    |                           |               |       vn[2]| |vn[3]       |
    |                           |               |            | |            |
    |                           |               |   son[0]   | |   son[1]   |
    |                           |               |            | |            |
    |                           |               |            | |            |
    |                           |               *------------* *------------*
    *-------------*-------------*
    vn[0]        vn[1] vn[0]        vn[1]
    vn[0]           en[0]           vn[1]

    refinement '1':                             refinement '2':

    vn[3]                           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

    *---------------------------*               *------------* *------------*
    |                           |               |            | |            |
    |                           |               |            | |            |
    |          son[1]           |               |            | |            |
    |                           |               |            | |            |
    |                           |               |            | |            |
    vn[0] *---------------------------* vn[1]         |            | |            |
    vn[3] *---------------------------* vn[2]         |   son[2]   | |   son[3]   |
    |                           |               |            | |            |
    |                           |               |            | |            |
    |          son[0]           |               |            | |            |
    |                           |               |            | |            |
    |                           |               |            | |            |
    *---------------------------*               *------------* *------------*

    vn[0]                           vn[1]       vn[0]        vn[1] vn[0]        vn[1]
    */

    /// Helper macros for easy iteration through all elements, nodes etc. in a Mesh.
#define for_all_elements(e, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
  if(((e) = (mesh)->get_element_fast(_id))->used)

#define for_all_base_elements(e, mesh) \
  for (int _id = 0; _id < (mesh)->get_num_base_elements(); _id++) \
  if(((e) = (mesh)->get_element_fast(_id))->used)

#define for_all_base_elements_incl_inactive(e, mesh) \
  for (int _id = 0; _id < (mesh)->get_num_base_elements(); _id++) \
  if(((e) = (mesh)->get_element_fast(_id))->used || !((e) = (mesh)->get_element_fast(_id))->used)

#define for_all_active_elements(e, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
  if(((e) = (mesh)->get_element_fast(_id))->used) \
  if((e)->active)

#define for_all_inactive_elements(e, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
  if(((e) = (mesh)->get_element_fast(_id))->used) \
  if(!(e)->active)

#define for_all_nodes(n, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
  if(((n) = (mesh)->get_node(_id))->used)

#define for_all_vertex_nodes(n, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
  if(((n) = (mesh)->get_node(_id))->used) \
  if(!(n)->type)

#define for_all_edge_nodes(n, mesh) \
  for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
  if(((n) = (mesh)->get_node(_id))->used) \
  if((n)->type)
  }
}
#endif
