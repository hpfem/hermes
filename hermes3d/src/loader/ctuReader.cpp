// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "ctuReader.h"
#include <string.h>
#include "../../../hermes_common/error.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/callstack.h"
#include "../mesh.h"
#include "../refdomain.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstdlib>

// maximal row length in bytes (used for reading the mesh3d-file)
#define MAX_ROW_LEN	                        1024

// exception error codes
#define E_CANT_OPEN_FILE                    -1
#define E_READ_ERROR                        -2

// number of markers on mesd3d file
#define MARKERS                             1

using namespace std;

struct _Node_
{
    double n[3];
};

struct _Hex_
{
    unsigned int n[8];
};

struct _Quad_
{
    unsigned int n[4];
};

class CTUInfo
{
    public:
        vector<_Node_*> nodes;
        vector<_Hex_*> hexs; 
        vector<_Quad_*> quads;

        CTUInfo()
        {
            
        }

        ~CTUInfo()
        {
            for(unsigned i=0;i<nodes.size();i++)
            {
                delete[] nodes.at(i);
            }
            nodes.clear();

            for(unsigned i=0;i<hexs.size();i++)
            {
                delete[] hexs.at(i);
            }
            hexs.clear();

            for(unsigned i=0;i<quads.size();i++)
            {
                delete[] quads.at(i);
            }
            quads.clear();
        }
};

bool in_list(vector<int> *l, int i)
{
    vector<int>::iterator it;

    for(it = l->begin(); it != l->end(); it++)
    {
        if(*it == i)
        {
            return true;
        }
    }

    return false;
}

void Split(std::vector<std::string>& lst, const std::string& input, const std::string& separators, bool remove_empty = true)
{
    string word = "";
    for (size_t n = 0; n < input.size(); ++n)
    {
        if (string::npos == separators.find(input[n]))
            word += input[n];
        else
        {
            if (!word.empty() || !remove_empty)
                lst.push_back(word);
            word = "";
        }
    }
    if (!word.empty() || !remove_empty)
        lst.push_back(word);
}

void Trim(string& str)
{
    // trim leading spaces
    size_t startpos = str.find_first_not_of(" \n\r\t");
    if( string::npos != startpos )
    {
        str = str.substr( startpos );
    }

    // trim trailing spaces
    size_t endpos = str.find_last_not_of(" \n\r\t");
    if( string::npos != endpos )
    {
        str = str.substr( 0, endpos+1 );
    }
}

void parse_ctuFormat(const char *file_name, CTUInfo *ctuInfo)
{
    ifstream ctu_f(file_name);

    bool point_start = true;
    bool point_write = true;

    string line, str;
    vector<string> cols;
    stringstream caster;

    while(getline(ctu_f,line))
    {
        Trim(line);

        Split(cols, line, " ",true);

        if(cols.size() == 1)
        {
            if(point_start)
            {
                point_start = false;
            }
            else
            {
                point_write = false;
            }
        }
        else
        {
            if(point_write)
            {
                _Node_ *node = new _Node_;

		node->n[0]  = atof(cols[1].c_str());
		node->n[1]  = atof(cols[2].c_str());
		node->n[2] = atof(cols[3].c_str());

                ctuInfo->nodes.push_back(node);
            }
            else
            {
                _Hex_ *hex = new _Hex_;

		hex->n[0]  = atof(cols[2].c_str());
		hex->n[1]  = atof(cols[3].c_str());
		hex->n[2]  = atof(cols[4].c_str());
		hex->n[3]  = atof(cols[5].c_str());
		hex->n[4]  = atof(cols[6].c_str());
		hex->n[5]  = atof(cols[7].c_str());
		hex->n[6]  = atof(cols[8].c_str());
		hex->n[7]  = atof(cols[9].c_str());

                ctuInfo->hexs.push_back(hex);
            }
        }

        cols.clear();
    }
    ctu_f.close();
}

CTUReader::CTUReader() {
    _F_
}

CTUReader::~CTUReader() {
    _F_
}

bool CTUReader::load(const char *file_name, Mesh *mesh)
{
    _F_
    assert(mesh != NULL);

    CTUInfo ci;
    parse_ctuFormat(file_name, &ci);

    vector<_Node_*>::iterator itv;
    vector<_Hex_*>::iterator ith;

    for(itv = ci.nodes.begin();itv != ci.nodes.end(); itv++)
    {
        mesh->add_vertex((*itv)->n[0], (*itv)->n[1], (*itv)->n[2]);
    }

    for(ith = ci.hexs.begin();ith < ci.hexs.end(); ith++)
    {
        mesh->add_hex((*ith)->n);
    }

    // check if all "outer" faces have defined boundary condition
    for (std::map<Facet::Key, Facet*>::const_iterator it = mesh->facets.begin(); it != mesh->facets.end(); it++)
    {
        Facet *facet = it->second;

        if(((unsigned) facet->left == INVALID_IDX) || ((unsigned) facet->right == INVALID_IDX))
        {
            fprintf(stderr, "Not all outer faces have defined boundary condition (line %d).", line_nr);
            throw E_READ_ERROR;
        }
    }

    mesh->ugh();
    return true;
}

bool CTUReader::save_as_h3d(const char *file_name, Mesh *mesh)
{
	_F_
	assert(mesh != NULL);

	FILE *file = fopen(file_name, "w");
	if (file == NULL) return false;

	// save vertices
	fprintf(file, "# vertices\n");
	fprintf(file, "%lu\n", (unsigned long int)mesh->vertices.size());
	for(std::map<unsigned int, Vertex*>::const_iterator it = mesh->vertices.begin(); it != mesh->vertices.end(); it++) {
    Vertex *v = it->second;
		fprintf(file, "%lf %lf %lf\n", v->x, v->y, v->z);
	}
	fprintf(file, "\n");

	// elements
	std::map<unsigned int, Element *> tet, hex, pri;
	for(std::map<unsigned int, Element*>::const_iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++) {
    Element *elem = it->second;
		if (elem->active) {
			switch (elem->get_mode()) {
      case HERMES_MODE_TET: tet[it->first] = elem; break;
				case HERMES_MODE_HEX: hex[it->first] = elem; break;
				case HERMES_MODE_PRISM: pri[it->first] = elem; break;
			}
		}
	}

	// save tetras
	fprintf(file, "# tetras\n");
	fprintf(file, "%lu\n", (unsigned long int)tet.size());
  for(std::map<unsigned int, Element*>::const_iterator it = tet.begin(); it != tet.end(); it++) {
		unsigned int vtcs[Tetra::NUM_VERTICES];
    it->second->get_vertices(vtcs);
		fprintf(file, "%u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3]);
	}
	fprintf(file, "\n");

	// save hexes
	fprintf(file, "# hexes\n");
	fprintf(file, "%lu\n", (unsigned long int)hex.size());
	for(std::map<unsigned int, Element*>::const_iterator it = hex.begin(); it != hex.end(); it++) {
		unsigned int vtcs[Hex::NUM_VERTICES];
    it->second->get_vertices(vtcs);
		fprintf(file, "%u %u %u %u %u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5], vtcs[6], vtcs[7]);
	}
	fprintf(file, "\n");

	// save prisms
	fprintf(file, "# prisms\n");
	fprintf(file, "%lu\n", (unsigned long int)pri.size());
  for(std::map<unsigned int, Element*>::const_iterator it = pri.begin(); it != pri.end(); it++) {
		unsigned int vtcs[Prism::NUM_VERTICES];
		it->second->get_vertices(vtcs);
		fprintf(file, "%u %u %u %u %u %u\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], vtcs[4], vtcs[5]);
	}
	fprintf(file, "\n");

	// boundaries
	std::map<unsigned int, Facet *> tri_facets, quad_facets;
  for(std::map<Facet::Key, Facet*>::iterator it = mesh->facets.begin(); it != mesh->facets.end(); it++) {
    Facet *facet = it->second;
		if(facet->type == Facet::OUTER && mesh->elements[facet->left]->active) {
			switch (facet->type) {
				case HERMES_MODE_TRIANGLE: 
          unsigned int ii;
          for(ii = 0; ; ii++)
            if(tri_facets[ii] == NULL)
              break;
          tri_facets[ii] = facet;
          break;
				case HERMES_MODE_QUAD: 
          unsigned int ij;
          for(ij = 0; ; ij++)
            if(quad_facets[ij] == NULL)
              break;
          quad_facets[ij] = facet;
          break;
			}
		}
	}

	// tris
	fprintf(file, "# tris\n");
	fprintf(file, "%lu\n", (unsigned long int)tri_facets.size());
	for(std::map<unsigned int, Facet*>::const_iterator it = tri_facets.begin(); it != tri_facets.end(); it++) {
    Facet *facet = it->second;
		Boundary *bnd = mesh->boundaries[facet->right];
		Element *elem = mesh->elements[facet->left];

		unsigned int vtcs[Tri::NUM_VERTICES];
		elem->get_face_vertices(facet->left_face_num, vtcs);

		fprintf(file, "%u %u %u     %d\n", vtcs[0], vtcs[1], vtcs[2], bnd->marker);
	}
	fprintf(file, "\n");

	// quads
	fprintf(file, "# quads\n");
	fprintf(file, "%lu\n", (unsigned long int)quad_facets.size());
	for(std::map<unsigned int, Facet*>::const_iterator it = quad_facets.begin(); it != quad_facets.end(); it++) {
    Facet *facet = it->second;
		Boundary *bnd = mesh->boundaries[facet->right];
		Element *elem = mesh->elements[facet->left];

		unsigned int vtcs[Quad::NUM_VERTICES];
		elem->get_face_vertices(facet->left_face_num, vtcs);

		fprintf(file, "%u %u %u %u     %d\n", vtcs[0], vtcs[1], vtcs[2], vtcs[3], bnd->marker);
	}
	fprintf(file, "\n");

	// the end
	fclose(file);

	return true;
}

bool CTUReader::save(const char *file_name, Mesh *mesh)
{
    return save_as_h3d(file_name, mesh);
}
