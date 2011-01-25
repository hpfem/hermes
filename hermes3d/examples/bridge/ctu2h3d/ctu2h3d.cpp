#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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

    while(getline(ctu_f,line))
    {
        boost::trim(line);

        boost::split(cols, line, boost::is_any_of("\t "));

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

                if(boost::lexical_cast<int>(cols[cols.size()-1]) >= 10 && boost::lexical_cast<int>(cols[cols.size()-1]) <= 16)
                {
                    //pd0_f << cols[0] << " " << cols[1] << " " << cols[2] << " " << cols[3] << "\n";
                    pd0_f << cols[1] << " " << cols[2] << " " << cols[3] << "\n";
                    pd0_id.push_back(boost::lexical_cast<int>(cols[0]));
                }
            }
            else
            {
                h3d_f << cols[2] << " " + cols[3] << " " << cols[4] << " " << cols[5] << " " << cols[6] << " " << cols[7] << " " << cols[8] << " " << cols[9] << "\n";

                str = cols[0] + " " + cols[1] + " ";
                bool hex_with_disp0 = false;

                for(int i=2;i<10;i++)
                {
                    if(in_list(&pd0_id, boost::lexical_cast<int>(cols[i])))
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
    }

    h3d_f << "\n# prims\n0\n\n# tri\n0\n\n";

    h3d_f.close();
    pd0_f.close();
    hd0_f.close();
    ctu_f.close();

return 0;
}
