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

    ctu_f.close();
}

void make_quad(unsigned n[], CTUInfo *ci)
{
    _Quad_ *ql[6];

    for(int i=0;i<6;i++)
        ql[i] = new _Quad_;

    ql[0]->n[0] = n[0];
    ql[0]->n[1] = n[1];
    ql[0]->n[2] = n[2];
    ql[0]->n[3] = n[3];

    ql[1]->n[0] = n[4];
    ql[1]->n[1] = n[5];
    ql[1]->n[2] = n[6];
    ql[1]->n[3] = n[7];

    ql[2]->n[0] = n[0];
    ql[2]->n[1] = n[1];
    ql[2]->n[2] = n[5];
    ql[2]->n[3] = n[4];

    ql[3]->n[0] = n[3];
    ql[3]->n[1] = n[2];
    ql[3]->n[2] = n[6];
    ql[3]->n[3] = n[7];

    ql[4]->n[0] = n[0];
    ql[4]->n[1] = n[3];
    ql[4]->n[2] = n[7];
    ql[4]->n[3] = n[4];

    ql[5]->n[0] = n[1];
    ql[5]->n[1] = n[2];
    ql[5]->n[2] = n[6];
    ql[5]->n[3] = n[5];

    for(int i=0;i<6;i++)
    {
        ci->quads.push_back(ql[i]);
    }
}

int main()
{
    CTUInfo ci;
    parse_ctuFormat("../most-sup.top", &ci);

    ofstream h3d_f("mesh.mesh3d");
    ofstream p_dbc_f("points_with_DBC.txt");
    ofstream h_dbc_f("hexs_with_DBC.txt");

    p_dbc_f << "id x y z" << endl;
    h_dbc_f << "id n1 n2 n3 n4 n5 n6 n7 n8" << endl;
    h3d_f << "# vertices" << endl << ci.nodes.size() << endl;

    cout << "Generating points for mesh.mesh3d & points_with_DBC.txt ..." << endl;

    for(int i=0;i<ci.nodes.size();i++)
    {
        _Node_ *n = ci.nodes.at(i);

        h3d_f << n->n[0] << " " << n->n[1] << " " << n->n[2] << endl;

        if(n->zero_DBC)
            p_dbc_f << n->id << " " << n->n[0] << " " << n->n[1] << " " << n->n[2] << endl;
    }

    h3d_f << endl << endl << "# tetras" << endl << 0 << endl;

    h3d_f << endl << "# hexes" << endl << ci.hexs.size() << endl;
    bool write_hex = true;

    cout << "Generating HEXS for mesh.mesh3d & hexs_with_DBC.txt ..." << endl;

    for(int i=0;i<ci.nodes.size();i++)
    {
        _Node_ *n = ci.nodes.at(i);

        for(int j=0;j<ci.hexs.size();j++)
        {
            _Hex_ *h = ci.hexs.at(j);

            for(int k=0;k<8;k++)
            {
                if(n->id == h->n[k] && n->zero_DBC == true)
                {
                    h->zero_DBC = true;
                    make_quad(h->n, &ci);

                    h_dbc_f << h->id << " " << h->n[0] << " " << h->n[1] << " " << h->n[2] << " " << h->n[3] << " " << h->n[4] << " " << h->n[5] << " " << h->n[6] << " " << h->n[7] << endl;
                    break;
                }
            }

            if(write_hex)
                h3d_f << h->n[0] << " " << h->n[1] << " " << h->n[2] << " " << h->n[3] << " " << h->n[4] << " " << h->n[5] << " " << h->n[6] << " " << h->n[7] << endl;
        }

        write_hex = false;
        //break;
    }

    h3d_f << endl << "# prisms" << endl << 0 << endl << endl;
    h3d_f << "# tris" << endl << 0 << endl << endl;

    cout << "Generating QUADS for mesh.mesh3d ..." << endl;
    h3d_f << endl << "# quads" << endl << ci.quads.size() << endl;

    for(int i=0;i<ci.quads.size();i++)
    {
        _Quad_ *q = ci.quads.at(i);
        h3d_f << q->n[0] << " " << q->n[1] << " " << q->n[2] << " " << q->n[3] << endl;
    }

    h3d_f.close();
    p_dbc_f.close();
    h_dbc_f.close();

return 0;
}
