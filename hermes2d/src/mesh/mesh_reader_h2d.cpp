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

#include "mesh.h"
#include "mesh_reader_h2d.h"
#include "refmap.h"

namespace Hermes
{
	namespace Hermes2D
	{
		extern unsigned g_mesh_seq;

		MeshReaderH2D::MeshReaderH2D()
		{
		}

		MeshReaderH2D::~MeshReaderH2D()
		{
		}

		Arc* MeshReaderH2D::load_arc(MeshSharedPtr mesh, MeshData *m, int id, Node** en, int &p1, int &p2, Arc* arc)
		{
			// read the end point indices
			p1 = m->curv_first[id];
			p2 = m->curv_second[id];

			*en = mesh->peek_edge_node(p1, p2);
			if (*en == nullptr)
				throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: edge %d-%d does not exist.", id, p1, p2);

			// get the number of control points
			std::vector<double> vCP;

			// edge endpoints are also control points, with weight 1.0
			arc->pt[0][0] = mesh->nodes[p1].x;
			arc->pt[0][1] = mesh->nodes[p1].y;
			arc->pt[0][2] = 1.0;
			arc->pt[2][0] = mesh->nodes[p2].x;
			arc->pt[2][1] = mesh->nodes[p2].y;
			arc->pt[2][2] = 1.0;

			// read the arc angle
			arc->angle = m->curv_third[id];
			double a = (180.0 - arc->angle) / 180.0 * M_PI;

			// generate one control point
			double x = 1.0 / tan(a * 0.5);
			arc->pt[1][0] = 0.5*((arc->pt[2][0] + arc->pt[0][0]) + (arc->pt[2][1] - arc->pt[0][1]) * x);
			arc->pt[1][1] = 0.5*((arc->pt[2][1] + arc->pt[0][1]) - (arc->pt[2][0] - arc->pt[0][0]) * x);
			arc->pt[1][2] = cos((M_PI - a) * 0.5);

			// get the number of knot vector points
			std::vector<double> vKnots;

			return arc;
		}

		Nurbs* MeshReaderH2D::load_nurbs(MeshSharedPtr mesh, MeshData *m, int id, Node** en, int &p1, int &p2, Nurbs* nurbs)
		{
			// read the end point indices
			p1 = m->curv_first[id];
			p2 = m->curv_second[id];

			*en = mesh->peek_edge_node(p1, p2);
			if (*en == nullptr)
				throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: edge %d-%d does not exist.", id, p1, p2);

			nurbs->degree = m->curv_third[id];

			// get the number of control points
			std::vector<double> vCP;

			int inner = 1, outer;
			for (unsigned int i = 0; i < m->vars_[m->curv_inner_pts[id]].size(); ++i)
			{
				std::istringstream istr(m->vars_[m->curv_inner_pts[id]][i]);

				double dummy_dbl;
				if (!(istr >> dummy_dbl))
					vCP.push_back(atof(m->vars_[m->vars_[m->curv_inner_pts[id]][i]][0].c_str()));
				else
					vCP.push_back(atof(m->vars_[m->curv_inner_pts[id]][i].c_str()));
			}
			inner = vCP.size() / 3;

			nurbs->np = inner + 2;

			// edge endpoints are also control points, with weight 1.0
			nurbs->pt = new double3[nurbs->np];
			nurbs->pt[0][0] = mesh->nodes[p1].x;
			nurbs->pt[0][1] = mesh->nodes[p1].y;
			nurbs->pt[0][2] = 1.0;
			nurbs->pt[inner + 1][0] = mesh->nodes[p2].x;
			nurbs->pt[inner + 1][1] = mesh->nodes[p2].y;
			nurbs->pt[inner + 1][2] = 1.0;

			// read inner control points
			for (int i = 0; i < inner; i++)
			{
				for (int j = 0; j < 3; ++j)
					nurbs->pt[i + 1][j] = vCP[3 * i + j];
			}

			// get the number of knot vector points
			std::vector<double> vKnots;

			inner = 0;
			for (unsigned int i = 0; i < m->vars_[m->curv_knots[id]].size(); ++i)
			{
				std::istringstream istr(m->vars_[m->curv_knots[id]][i]);
				double dummy_dbl;

				if (!(istr >> dummy_dbl))
					vKnots.push_back(atof(m->vars_[m->vars_[m->curv_knots[id]][i]][0].c_str()));
				else
					vKnots.push_back(atof(m->vars_[m->curv_knots[id]][i].c_str()));
			}
			inner = vKnots.size();

			nurbs->nk = nurbs->degree + nurbs->np + 1;
			outer = nurbs->nk - inner;
			if ((outer & 1) == 1)
				throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: incorrect number of knot points.", id);

			// knot vector is completed by 0.0 on the left and by 1.0 on the right
			nurbs->kv = new double[nurbs->nk];

			for (int i = 0; i < outer / 2; i++)
				nurbs->kv[i] = 0.0;

			if (inner) {
				for (int i = outer / 2; i < inner + outer / 2; i++) {
					nurbs->kv[i] = vKnots[i - (outer / 2)];
				}
			}

			for (int i = outer / 2 + inner; i < nurbs->nk; i++)
				nurbs->kv[i] = 1.0;

			return nurbs;
		}

		Curve* MeshReaderH2D::load_curve(MeshSharedPtr mesh, MeshData *m, int id, Node** en, int &p1, int &p2)
		{
			if (m->curv_nurbs[id])
				return load_nurbs(mesh, m, id, en, p1, p2, new Nurbs);
			else
				return load_arc(mesh, m, id, en, p1, p2, new Arc);
		}

		void MeshReaderH2D::load(std::string filename, MeshSharedPtr mesh)
		{
			// Check if file exists
			std::ifstream s(filename);
			if (!s.good())
				throw Hermes::Exceptions::MeshLoadFailureException("Mesh file not found.");
			s.close();

			int i, j, n;
			Node* en;
			bool debug = false;

			mesh->free();

			MeshData m(filename);
			m.parse_mesh();

			//// vertices ////////////////////////////////////////////////////////////////

			n = m.n_vert;
			if (n < 0) throw Hermes::Exceptions::MeshLoadFailureException("File %s: 'vertices' must be a list.", filename.c_str());
			if (n < 2) throw Hermes::Exceptions::MeshLoadFailureException("File %s: invalid number of vertices.", filename.c_str());

			// create a hash table large enough
			int size = HashTable::H2D_DEFAULT_HASH_SIZE;
			while (size < 8 * n) size *= 2;
			mesh->init(size);

			// create top-level vertex nodes
			for (i = 0; i < n; i++)
			{
				Node* node = mesh->nodes.add();
				assert(node->id == i);
				node->ref = TOP_LEVEL_REF;
				node->type = HERMES_TYPE_VERTEX;
				node->bnd = 0;
				node->p1 = node->p2 = -1;
				node->next_hash = nullptr;
				node->x = m.x_vertex[i];
				node->y = m.y_vertex[i];
			}
			mesh->ntopvert = n;

			//// elements ////////////////////////////////////////////////////////////////

			n = m.n_el;
			if (n < 0) throw Hermes::Exceptions::MeshLoadFailureException("File %s: 'elements' must be a list.", filename.c_str());
			if (n < 1) throw Hermes::Exceptions::MeshLoadFailureException("File %s: no elements defined.", filename.c_str());

			// create elements
			mesh->nactive = 0;
			for (i = 0; i < n; i++)
			{
				// read and check vertex indices
				int nv;
				if (m.en4[i] == -1)  nv = 4;
				else nv = 5;

				int* idx = new int[nv - 1];
				std::string el_marker;
				if (!nv) {
					mesh->elements.skip_slot()->cm = nullptr;
					continue;
				}

				if (nv < 4 || nv > 5)
				{
					delete[] idx;
					throw Hermes::Exceptions::MeshLoadFailureException("File %s: element #%d: wrong number of vertex indices.", filename.c_str(), i);
				}

				if (nv == 4) {
					idx[0] = m.en1[i];
					idx[1] = m.en2[i];
					idx[2] = m.en3[i];

					el_marker = m.e_mtl[i];
				}
				else {
					idx[0] = m.en1[i];
					idx[1] = m.en2[i];
					idx[2] = m.en3[i];
					idx[3] = m.en4[i];

					el_marker = m.e_mtl[i];
				}
				for (j = 0; j < nv - 1; j++)
					if (idx[j] < 0 || idx[j] >= mesh->ntopvert)
					{
						delete[] idx;
						throw Hermes::Exceptions::MeshLoadFailureException("File %s: error creating element #%d: vertex #%d does not exist.", filename.c_str(), i, idx[j]);
					}

				Node *v0 = &mesh->nodes[idx[0]], *v1 = &mesh->nodes[idx[1]], *v2 = &mesh->nodes[idx[2]];

				int marker;

				// This functions check if the user-supplied marker on this element has been
				// already used, and if not, inserts it in the appropriate structure.
				mesh->element_markers_conversion.insert_marker(el_marker);
				marker = mesh->element_markers_conversion.get_internal_marker(el_marker).marker;

				if (nv == 4) {
					Mesh::check_triangle(i, v0, v1, v2);
					mesh->create_triangle(marker, v0, v1, v2, nullptr);
				}
				else {
					Node *v3 = &mesh->nodes[idx[3]];
					Mesh::check_quad(i, v0, v1, v2, v3);
					mesh->create_quad(marker, v0, v1, v2, v3, nullptr);
				}

				mesh->nactive++;

				delete[] idx;
			}
			mesh->nbase = n;

			//// boundaries //////////////////////////////////////////////////////////////
			if (m.n_bdy > 0)
			{
				n = m.n_bdy;

				// read boundary data
				for (i = 0; i < n; i++)
				{
					int v1, v2, marker;
					v1 = m.bdy_first[i];
					v2 = m.bdy_second[i];

					en = mesh->peek_edge_node(v1, v2);
					if (en == nullptr)
						throw Hermes::Exceptions::MeshLoadFailureException("File %s: boundary data #%d: edge %d-%d does not exist", filename.c_str(), i, v1, v2);

					std::string bnd_marker;
					bnd_marker = m.bdy_type[i];

					// This functions check if the user-supplied marker on this element has been
					// already used, and if not, inserts it in the appropriate structure.
					mesh->boundary_markers_conversion.insert_marker(bnd_marker);
					marker = mesh->boundary_markers_conversion.get_internal_marker(bnd_marker).marker;

					en->marker = marker;
				}
			}

			Node* node;
			for_all_edge_nodes(node, mesh)
			{
				if (node->ref < 2)
				{
					mesh->nodes[node->p1].bnd = 1;
					mesh->nodes[node->p2].bnd = 1;
					node->bnd = 1;
				}
			}

			//// curves //////////////////////////////////////////////////////////////////
			if (m.n_curv > 0)
			{
				n = m.n_curv;
				if (n < 0) throw Hermes::Exceptions::MeshLoadFailureException("File %s: 'curves' must be a list.", filename.c_str());

				// load curved edges
				for (i = 0; i < n; i++)
				{
					// load the control points, knot vector, etc.
					Node* en;
					int p1, p2;

					Curve* curve = load_curve(mesh, &m, i, &en, p1, p2);

					// assign the arc to the elements sharing the edge node
					MeshUtil::assign_curve(en, curve, p1, p2);
				}
			}

			// update refmap coeffs of curvilinear elements
			Element* e;
			for_all_used_elements(e, mesh)
			{
				if (e->cm != nullptr)
					e->cm->update_refmap_coeffs(e);
				this->ref_map.set_element_iro_cache(e);
			}

			//// refinements /////////////////////////////////////////////////////////////
			if (m.n_ref > 0)
			{
				n = m.n_ref;
				if (n < 0) throw Hermes::Exceptions::MeshLoadFailureException("File %s: 'refinements' must be a list.", filename.c_str());

				// perform initial refinements
				for (i = 0; i < n; i++)
				{
					int id, ref;
					id = m.ref_elt[i];
					ref = m.ref_type[i];
					mesh->refine_element_id(id, ref);
				}
			}
			mesh->ninitial = mesh->elements.get_num_items();

			mesh->seq = g_mesh_seq++;
			if (HermesCommonApi.get_integral_param_value(checkMeshesOnLoad))
				mesh->initial_single_check();
		}

		void MeshReaderH2D::load(const char *filename, MeshSharedPtr mesh)
		{
			std::string s(filename);
			return this->load(s, mesh);
		}

		void MeshReaderH2D::save_refinements(MeshSharedPtr mesh, FILE* f, Element* e, int id, bool& first)
		{
			if (e->active) return;
			fprintf(f, first ? "refinements =\n{\n" : ",\n"); first = false;
			if (e->bsplit())
			{
				fprintf(f, "  { %d, 0 }", id);
				int sid = mesh->seq; mesh->seq += 4;
				for (int i = 0; i < 4; i++)
					save_refinements(mesh, f, e->sons[i], sid + i, first);
			}
			else if (e->hsplit())
			{
				fprintf(f, "  { %d, 1 }", id);
				int sid = mesh->seq; mesh->seq += 2;
				save_refinements(mesh, f, e->sons[0], sid, first);
				save_refinements(mesh, f, e->sons[1], sid + 1, first);
			}
			else
			{
				fprintf(f, "  { %d, 2 }", id);
				int sid = mesh->seq; mesh->seq += 2;
				save_refinements(mesh, f, e->sons[2], sid, first);
				save_refinements(mesh, f, e->sons[3], sid + 1, first);
			}
		}

		void MeshReaderH2D::save_curve(MeshSharedPtr mesh, FILE* f, int p1, int p2, Curve* curve)
		{
			if (curve->type == ArcType)
			{
				fprintf(f, " [ %d, %d, %.16g ]", p1, p2, ((Arc*)curve)->angle);
			}
			else
			{
				int inner = ((Nurbs*)curve)->np - 2;
				int outer = ((Nurbs*)curve)->nk - inner;
				fprintf(f, "  ú %d, %d, %d, ú ", p1, p2, ((Nurbs*)curve)->degree);
				for (int i = 1; i < ((Nurbs*)curve)->np - 1; i++)
					fprintf(f, "ú %.16g, %.16g, %.16g ]%s ",
						((Nurbs*)curve)->pt[i][0], ((Nurbs*)curve)->pt[i][1], ((Nurbs*)curve)->pt[i][2],
						i < ((Nurbs*)curve)->np - 2 ? "," : "");

				fprintf(f, "],[ ");
				int max = ((Nurbs*)curve)->nk - (((Nurbs*)curve)->degree + 1);
				for (int i = ((Nurbs*)curve)->degree + 1; i < max; i++)
					fprintf(f, "%.16g%s", ((Nurbs*)curve)->kv[i], i < max - 1 ? "," : "");
				fprintf(f, "] ]");
			}
		}

		void MeshReaderH2D::save(const char* filename, MeshSharedPtr mesh)
		{
			std::string s(filename);
			return this->save(s, mesh);
		}

		void MeshReaderH2D::save(std::string filename, MeshSharedPtr mesh)
		{
			int i, mrk;
			Element* e;

			// open output file
			FILE* f = fopen(filename.c_str(), "w");
			if (f == nullptr) throw Hermes::Exceptions::MeshLoadFailureException("Could not create mesh file.");
			//fprintf(f, "# hermes2d saved mesh\n\n");

			// save vertices
			fprintf(f, "vertices =\n[\n");
			for (i = 0; i < mesh->ntopvert; i++)
				fprintf(f, " [ %.16g, %.16g ]%s\n", mesh->nodes[i].x, mesh->nodes[i].y, (i < mesh->ntopvert - 1 ? "," : ""));

			// save elements
			fprintf(f, "]\n\nelements =\n[");
			bool first = true;
			for (i = 0; i < mesh->get_num_base_elements(); i++)
			{
				const char* nl = first ? "\n" : ",\n";  first = false;
				e = mesh->get_element_fast(i);
				if (!e->used)
					fprintf(f, "%s [ ]", nl);
				else if (e->is_triangle())
					fprintf(f, "%s [ %d, %d, %d, \"%s\" ]", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str());
				else
					fprintf(f, "%s [ %d, %d, %d, %d, \"%s\" ]", nl, e->vn[0]->id, e->vn[1]->id, e->vn[2]->id, e->vn[3]->id, mesh->get_element_markers_conversion().get_user_marker(e->marker).marker.c_str());
			}

			// save boundary markers
			fprintf(f, "\n]\n\nboundaries =\n[");
			first = true;
			for_all_base_elements(e, mesh)
			{
				for (unsigned char i = 0; i < e->get_nvert(); i++)
					if ((mrk = MeshUtil::get_base_edge_node(e, i)->marker))
					{
						const char* nl = first ? "\n" : ",\n";  first = false;
						fprintf(f, "%s [ %d, %d, \"%s\" ]", nl, e->vn[i]->id, e->vn[e->next_vert(i)]->id, mesh->boundary_markers_conversion.get_user_marker(mrk).marker.c_str());
					}
			}
			fprintf(f, "\n]\n\n");

			// save curved edges
			first = true;
			for_all_base_elements(e, mesh)
			{
				if (e->is_curved())
				{
					for (unsigned char i = 0; i < e->get_nvert(); i++)
						if (e->cm->curves[i] != nullptr)
						{
							fprintf(f, first ? "curves =\n[\n" : ",\n");  first = false;
							save_curve(mesh, f, e->vn[i]->id, e->vn[e->next_vert(i)]->id, e->cm->curves[i]);
						}
					if (!first) fprintf(f, "\n]\n\n");
				}
			}
			// save refinements
			unsigned temp = mesh->seq;
			mesh->seq = mesh->nbase;
			first = true;
			for_all_base_elements(e, mesh)
				save_refinements(mesh, f, e, e->id, first);
			if (!first) fprintf(f, "\n]\n\n");

			mesh->seq = temp;
			fclose(f);
		}
	}
}