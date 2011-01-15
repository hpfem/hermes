#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

int main()
{
    ofstream o_f("op.txt");
    ifstream b_f("bridge.txt");
    ifstream s_f("supports.txt");

    string node_id, b_line;
    vector<string> cols;

        getline(b_f,b_line);
        getline(b_f,b_line);

        boost::trim(b_line);
        boost::split(cols, b_line, boost::is_any_of("\t "));

        cout << "col Size: " << cols.size() << endl;
        cout << cols[0] << " | " << cols[1] << " | " << cols[2] << " | " << cols[3] << " | " << cols[4] << " | " << cols[5] << " | " << endl;
/*
    while(getline(s_f,node_id))
    {
        boost::trim(node_id);
        cout << "Node id: " << node_id << endl;
        
        while(getline(b_f,b_line))
        {
            boost::trim(b_line);
            boost::split(cols, b_line, boost::is_any_of("\t  "));

            cout << "col Size: " << cols.size() << endl;

            break;

            if(cols.at(0) == node_id)
            {
                o_f << cols[1] << " " << cols[2] << " " << cols[3] << endl;
                break;
            }

            cols.clear();
        }

        if(b_f.eof())
            b_f.seekg(0,ios::beg);
    }
*/

    o_f.close();
    s_f.close();
    b_f.close();

    return 0;
}
