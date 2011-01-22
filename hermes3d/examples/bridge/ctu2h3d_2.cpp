#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

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

int main()
{
    ofstream h3d_f("test.mesh3d");
    ofstream pd0_f("test.pd0.mesh3d");
    ofstream hd0_f("test.hd0.mesh3d");
    ifstream ctu_f("most-sup.top");

    bool point_start = true;
    bool point_write = true;

    string line, str;

    vector<string> cols;
    vector<int> pd0_id;

    stringstream caster;

    while(getline(ctu_f,line))
    {
        Trim(line);

        Split(cols, line, " ",true);

        if(cols.size() == 1)
        {
            if(point_start)
            {
                h3d_f << "# points\n";
                pd0_f << "# points with displacement 0\n";
                point_start = false;
            }
            else
            {
                h3d_f << "# hex\n";
                hd0_f << "# hex with atleast one point with displacement 0\n";
                point_write = false;
            }

            h3d_f << cols[0] << "\n";
        }
        else
        {
            if(point_write)
            {
                h3d_f << cols[1] << " " << cols[2] << " " << cols[3] << "\n";

                int v;
                caster << cols[cols.size()-1];
                caster >> v;
                caster.str("");
                if(v >= 10 && v <= 16)
                {
                    pd0_f << cols[0] << " " << cols[1] << " " << cols[2] << " " << cols[3] << "\n";

                    int k;
                    caster << cols[0];
                    caster >> k;
                    caster.str("");

                    pd0_id.push_back(k);
                }
            }
            else
            {
                h3d_f << cols[2] << " " + cols[3] << " " << cols[4] << " " << cols[5] << " " << cols[6] << " " << cols[7] << " " << cols[8] << " " << cols[9] << "\n";

                str = cols[0] + " " + cols[1] + " ";
                bool hex_with_disp0 = false;

                for(int i=2;i<10;i++)
                {
                    int v1;
                    caster << cols[i];
                    caster >> v1;
                    caster.str("");
                    if(in_list(&pd0_id, v1))
                    {
                        str += "|" + cols[i] + "|" + " ";
                        hex_with_disp0 = true;
                    }
                    else
                        str += cols[i] + " ";
                }

                if(hex_with_disp0)
                    hd0_f << str + "\n";
            }
        }

        cols.clear();
    }

    h3d_f << "\n# prims\n0\n\n# tri\n0\n\n";

    h3d_f.close();
    pd0_f.close();
    hd0_f.close();
    ctu_f.close();

return 0;
}
