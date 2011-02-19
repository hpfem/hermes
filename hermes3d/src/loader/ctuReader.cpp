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
    unsigned id;
    double n[3];
    bool zero_DBC;
};

struct _Hex_
{
    unsigned id;
    unsigned int n[8];
    bool zero_DBC;
};

struct _Quad_
{
    unsigned id;
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
                delete nodes.at(i);
            }
            nodes.clear();

            for(unsigned i=0;i<hexs.size();i++)
            {
                delete hexs.at(i);
            }
            hexs.clear();

            for(unsigned i=0;i<quads.size();i++)
            {
                delete quads.at(i);
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

void make_quad(_Hex_ *hex, CTUInfo *ci)
{
    _Node_ *n1 = ci->nodes.at(hex->n[6]-1);
    _Node_ *n2 = ci->nodes.at(hex->n[7]-1);
    _Node_ *n3 = ci->nodes.at(hex->n[4]-1);
    _Node_ *n4 = ci->nodes.at(hex->n[5]-1);
    _Node_ *n5 = ci->nodes.at(hex->n[2]-1);
    _Node_ *n6 = ci->nodes.at(hex->n[3]-1);
    _Node_ *n7 = ci->nodes.at(hex->n[0]-1);
    _Node_ *n8 = ci->nodes.at(hex->n[1]-1);

    if (n1->zero_DBC == true && n2->zero_DBC == true && n3->zero_DBC == true && n4->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;

        quad->n[0] = hex->n[6];
        quad->n[1] = hex->n[7];
        quad->n[2] = hex->n[4];
        quad->n[3] = hex->n[5];
        
        ci->quads.push_back(quad);
    }

    if (n1->zero_DBC == true && n2->zero_DBC == true && n6->zero_DBC == true && n5->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;
        
        quad->n[0] = hex->n[6];
        quad->n[1] = hex->n[7];
        quad->n[2] = hex->n[3];
        quad->n[3] = hex->n[2];

        ci->quads.push_back(quad);
    }
    
    if (n2->zero_DBC == true && n3->zero_DBC == true && n7->zero_DBC == true && n6->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;
        
        quad->n[0] = hex->n[7];
        quad->n[1] = hex->n[4];
        quad->n[2] = hex->n[0];
        quad->n[3] = hex->n[3];

        ci->quads.push_back(quad);
    }

    if (n3->zero_DBC == true && n4->zero_DBC == true && n8->zero_DBC == true && n7->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;
        
        quad->n[0] = hex->n[4];
        quad->n[1] = hex->n[5];
        quad->n[2] = hex->n[1];
        quad->n[3] = hex->n[0];

        ci->quads.push_back(quad);
    }

    if (n4->zero_DBC == true && n1->zero_DBC == true && n5->zero_DBC == true && n8->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;
        
        quad->n[0] = hex->n[5];
        quad->n[1] = hex->n[6];
        quad->n[2] = hex->n[2];
        quad->n[3] = hex->n[1];

        ci->quads.push_back(quad);
    }

    if (n5->zero_DBC == true && n6->zero_DBC == true && n7->zero_DBC == true && n8->zero_DBC == true)
    {
        _Quad_ *quad = new _Quad_;
        
        quad->id = ci->quads.size()+1;
        
        quad->n[0] = hex->n[2];
        quad->n[1] = hex->n[3];
        quad->n[2] = hex->n[0];
        quad->n[3] = hex->n[1];

        ci->quads.push_back(quad);
    }
}

void parse_ctuFormat(const char *file_name, CTUInfo *ctuInfo)
{
    ifstream ctu_f(file_name);
    string line, str;
    vector<string> cols;

    //read number of vertices
    getline(ctu_f,line);

    Trim(line);
    int num_of_vertices = atoi(line.c_str());
    
    //read vertices
    for(int i=0;i<num_of_vertices;i++)
    {
        getline(ctu_f,line);
        Trim(line);
        Split(cols, line, " ",true);

        _Node_ *node = new _Node_;

        node->id = atoi(cols[0].c_str());
        node->n[0]  = atof(cols[1].c_str());
        node->n[1]  = atof(cols[2].c_str());
        node->n[2] = atof(cols[3].c_str());

        int c = atoi(cols[cols.size()-1].c_str());

        if( c >= 10 && c <= 16)
            node->zero_DBC = true;
        else
            node->zero_DBC = false;

        ctuInfo->nodes.push_back(node);
        cols.clear();
    }

    //read number of hexs
    getline(ctu_f,line);

    Trim(line);
    int num_of_hexs = atoi(line.c_str());

    //read vertices
    for(int i=0;i<num_of_hexs;i++)
    {
        getline(ctu_f,line);
        Trim(line);
        Split(cols, line, " ",true);

        _Hex_ *hex = new _Hex_;

        hex->id  = atoi(cols[0].c_str());
        hex->n[0]  = atof(cols[2].c_str());
        hex->n[1]  = atof(cols[3].c_str());
        hex->n[2]  = atof(cols[4].c_str());
        hex->n[3]  = atof(cols[5].c_str());
        hex->n[4]  = atof(cols[6].c_str());
        hex->n[5]  = atof(cols[7].c_str());
        hex->n[6]  = atof(cols[8].c_str());
        hex->n[7]  = atof(cols[9].c_str());
        hex->zero_DBC  = false;

        ctuInfo->hexs.push_back(hex);

        make_quad(hex, ctuInfo);

        cols.clear();
    }
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
    vector<_Quad_*>::iterator itq;

    for(itv = ci.nodes.begin();itv != ci.nodes.end(); itv++)
    {
        mesh->add_vertex((*itv)->n[0], (*itv)->n[1], (*itv)->n[2]);
    }

    for(ith = ci.hexs.begin();ith < ci.hexs.end(); ith++)
    {
        mesh->add_hex((*ith)->n);
    }

    for(itq = ci.quads.begin();itq < ci.quads.end(); itq++)
    {
		mesh->add_quad_boundary((*ith)->n, 1);
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
