#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

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

/*
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
            }
            else
            {
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
            }
        }

        cols.clear();
    }
}
*/

int main()
{
    CTUInfo ci;
    parse_ctuFormat("../data/most-sup.top", &ci);

    ofstream h3d_f("../data/most.mesh3d");

    h3d_f << "# vertices" << endl << ci.nodes.size() << endl;

    cout << "Generating VERTICES ..." << endl;

    for(unsigned int i=0;i<ci.nodes.size();i++)
    {
        _Node_ *n = ci.nodes.at(i);
        h3d_f << n->n[0] << " " << n->n[1] << " " << n->n[2] << endl;
    }

    h3d_f << endl << endl << "# tetras" << endl << 0 << endl;
    h3d_f << endl << "# hexes" << endl << ci.hexs.size() << endl;
    bool write_hex = true;

    cout << "Generating HEXS ..." << endl;

    for (unsigned int i = 0; i < ci.hexs.size(); i++)
    {
        _Hex_ *h = ci.hexs.at(i);
        h3d_f <<  h->n[6] << " " << h->n[7] << " " << h->n[4] << " " << h->n[5] << " " << h->n[2] << " " << h->n[3] << " " << h->n[0] << " " << h->n[1] << endl;
    }     

    h3d_f << endl << "# prisms" << endl << 0 << endl << endl;
    h3d_f << "# tris" << endl << 0 << endl;

    cout << "Generating QUADS ..." << endl;
    h3d_f << endl << "# quads" << endl << ci.quads.size() << endl;

    for (unsigned int i=0; i < ci.quads.size(); i++)
    {
        _Quad_ *q = ci.quads.at(i);
        h3d_f <<  q->n[0] << " " << q->n[1] << " " << q->n[2] << " " << q->n[3] << " 1" << endl;
    }

/*
    int nquads=0;
    for (int i=0;i<ci.hexs.size();i++)
    {
        _Hex_ *h = ci.hexs.at(i);
        _Node_ *n1 = ci.nodes.at(h->n[6]-1);
        _Node_ *n2 = ci.nodes.at(h->n[7]-1);
        _Node_ *n3 = ci.nodes.at(h->n[4]-1);
        _Node_ *n4 = ci.nodes.at(h->n[5]-1);
        _Node_ *n5 = ci.nodes.at(h->n[2]-1);
        _Node_ *n6 = ci.nodes.at(h->n[3]-1);
        _Node_ *n7 = ci.nodes.at(h->n[0]-1);
        _Node_ *n8 = ci.nodes.at(h->n[1]-1);
       
        if (n1->zero_DBC == true && n2->zero_DBC == true && n3->zero_DBC == true && n4->zero_DBC == true)
            nquads++;
        if (n1->zero_DBC == true && n2->zero_DBC == true && n6->zero_DBC == true && n5->zero_DBC == true)
            nquads++;
        if (n2->zero_DBC == true && n3->zero_DBC == true && n7->zero_DBC == true && n6->zero_DBC == true)
            nquads++;
        if (n3->zero_DBC == true && n4->zero_DBC == true && n8->zero_DBC == true && n7->zero_DBC == true)
            nquads++;
        if (n4->zero_DBC == true && n1->zero_DBC == true && n5->zero_DBC == true && n8->zero_DBC == true)
            nquads++;
        if (n5->zero_DBC == true && n6->zero_DBC == true && n7->zero_DBC == true && n8->zero_DBC == true)
            nquads++;
    }     

    for (int i=0;i<ci.hexs.size();i++)
    {
        _Hex_ *h = ci.hexs.at(i);
        _Node_ *n1 = ci.nodes.at(h->n[6]-1);
        _Node_ *n2 = ci.nodes.at(h->n[7]-1);
        _Node_ *n3 = ci.nodes.at(h->n[4]-1);
        _Node_ *n4 = ci.nodes.at(h->n[5]-1);
        _Node_ *n5 = ci.nodes.at(h->n[2]-1);
        _Node_ *n6 = ci.nodes.at(h->n[3]-1);
        _Node_ *n7 = ci.nodes.at(h->n[0]-1);
        _Node_ *n8 = ci.nodes.at(h->n[1]-1);

        if (n1->zero_DBC == true && n2->zero_DBC == true && n3->zero_DBC == true && n4->zero_DBC == true)
            h3d_f << h->n[6] << " " << h->n[7] << " " << h->n[4] << " " << h->n[5] << " 1" << endl;
        if (n1->zero_DBC == true && n2->zero_DBC == true && n6->zero_DBC == true && n5->zero_DBC == true)
            h3d_f << h->n[6] << " " << h->n[7] << " " << h->n[3] << " " << h->n[2] << " 1" << endl;
        if (n2->zero_DBC == true && n3->zero_DBC == true && n7->zero_DBC == true && n6->zero_DBC == true)
            h3d_f << h->n[7] << " " << h->n[4] << " " << h->n[0] << " " << h->n[3] << " 1" << endl;
        if (n3->zero_DBC == true && n4->zero_DBC == true && n8->zero_DBC == true && n7->zero_DBC == true)
            h3d_f << h->n[4] << " " << h->n[5] << " " << h->n[1] << " " << h->n[0] << " 1" << endl;
        if (n4->zero_DBC == true && n1->zero_DBC == true && n5->zero_DBC == true && n8->zero_DBC == true)
            h3d_f << h->n[5] << " " << h->n[6] << " " << h->n[2] << " " << h->n[1] << " 1" << endl;
        if (n5->zero_DBC == true && n6->zero_DBC == true && n7->zero_DBC == true && n8->zero_DBC == true)
            h3d_f << h->n[2] << " " << h->n[3] << " " << h->n[0] << " " << h->n[1] << " 1" << endl;
    }  
*/
return 0;
}
