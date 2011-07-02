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

# include "mesh2d_parser.h"
namespace Hermes
{
  namespace Hermes2D
  {
    void MeshData::strip(std::string& str)
    {
      std::string temp;

      // Remove comments
      if (str.find('#') != str.npos)
        str.erase(str.find('#'));

      // Remove brackets, commas and unnecessary blank spaces
      for (size_t i = 0; i < str.length(); i++)
      {
        if (str[i] != ' ' && str[i] != '\t' && str[i] != '[' && str[i] != ']' && str[i] != '{' && str[i] != '}' && str[i] != '"')
        {
          if (str[i] == ',' || str[i] == ';')
            temp.append("\t");
          else if (str[i] == '=')
            temp.append("=\t");
          else	
            temp.append(1,str[i]);
        }
      }

      str.assign(temp);
    }

    void MeshData::parse_mesh(void)
    {
      std::vector<std::string> varlist;

      int dummy_int;
      double dummy_dbl;
      std::string dummy_str;

      std::ifstream inFile(mesh_file_.c_str());
      std::string line, word, temp_word;

      int counter(0);
      bool isVert(false), isElt(false), isBdy(false), isCurv(false), isRef(false), isVar(false);

      while (std::getline(inFile,line))
      {		
        // Remove all comments, unnecessary blank spaces, commas and paranthesis
        strip(line);

        // Display stripped file
        //if (line.find_first_not_of("\t ") != line.npos)
        //	std::cout << line << std::endl;

        if (line.find_first_not_of("\t ") != line.npos)
        {			
          std::istringstream stream(line);
          stream >> word;

          if (*word.rbegin() == '=')
          {
            word.erase(word.size() - 1);

            if (word == "vertices")
            {
              isVert = true; 
              isElt = false; isBdy = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if (word == "elements")
            {
              isElt = true; 
              isVert = false; isBdy = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if (word == "boundaries")
            {
              isBdy = true; 
              isVert = false; isElt = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if (word == "curves")
            {
              isCurv = true; 
              isVert = false; isElt = false; isBdy = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if (word == "refinements")
            {
              isRef = true; 
              isVert = false; isElt = false; isBdy = false; isCurv = false; isVar = false;
              counter = -1;
            }
            else
            {
              isVar = true;
              isVert = false; isElt = false; isBdy = false; isCurv = false; isRef = false;
              counter = -1;

              temp_word = word;
            }
          }

          if (counter == -1)
            counter = 0;
          else
          {
            if (isVert)
            {
              std::istringstream istr(word);

              if (!(istr >> dummy_dbl))
                x_vertex.push_back(atof(vars_[word][0].c_str()));
              else
                x_vertex.push_back(atof(word.c_str()));

              ++counter;
            }
            if (isElt)
            {
              std::istringstream istr(word);

              if (!(istr >> dummy_int))
                en1.push_back(atoi(vars_[word][0].c_str()));
              else
                en1.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if (isBdy)
            {
              std::istringstream istr(word);

              if (!(istr >> dummy_dbl))
                bdy_first.push_back(atof(vars_[word][0].c_str()));
              else
                bdy_first.push_back(atof(word.c_str()));					

              ++counter;
            }
            else if (isCurv)
            {
              std::istringstream istr(word);

              if (!(istr >> dummy_int))
                curv_first.push_back(atoi(vars_[word][0].c_str()));
              else
                curv_first.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if (isRef)
            {
              std::istringstream istr(word);

              if (!(istr >> dummy_int))
                ref_elt.push_back(atoi(vars_[word][0].c_str()));
              else 
                ref_elt.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if (isVar)
            {
              vars_[temp_word].push_back(word);
              ++counter;
            }
          }	

          if (isVert)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if (counter%2 == 0)
              {	
                if (!(istr >> dummy_dbl))
                  x_vertex.push_back(atof(vars_[word][0].c_str()));
                else
                  x_vertex.push_back(atof(word.c_str()));
              }
              else
              {	
                if (!(istr >> dummy_dbl))
                  y_vertex.push_back(atof(vars_[word][0].c_str()));
                else
                  y_vertex.push_back(atof(word.c_str()));
              }

              ++counter;
            }
          }
          else if (isElt)		
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if (counter%5 == 0)
              {
                if (!(istr >> dummy_int))
                  en1.push_back(atoi(vars_[word][0].c_str()));
                else
                  en1.push_back(atoi(word.c_str()));
              }
              else if (counter%5 == 1)
              {
                if (!(istr >> dummy_int))
                  en2.push_back(atoi(vars_[word][0].c_str()));
                else
                  en2.push_back(atoi(word.c_str()));
              }
              else if (counter%5 == 2)
              {
                if (!(istr >> dummy_int))
                  en3.push_back(atoi(vars_[word][0].c_str()));
                else
                  en3.push_back(atoi(word.c_str()));
              }
              else if (counter%5 == 3)
              {
                if (!(istr >> dummy_int))
                {	
                  en4.push_back(-1);
                  e_mtl.push_back(word);

                  ++counter;
                }
                else
                  en4.push_back(atoi(word.c_str()));
              }
              else
                e_mtl.push_back(word);

              ++counter;
            }
          }
          else if (isBdy)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if (counter%3 == 0)
              {
                if (!(istr >> dummy_int))
                  bdy_first.push_back(atoi(vars_[word][0].c_str()));
                else
                  bdy_first.push_back(atoi(word.c_str()));
              }
              else if (counter%3 == 1)
              {
                if (!(istr >> dummy_int))
                  bdy_second.push_back(atoi(vars_[word][0].c_str()));
                else
                  bdy_second.push_back(atoi(word.c_str()));
              }
              else
                bdy_type.push_back(word);

              ++counter;
            }
          }
          else if (isCurv)
          {
            while (stream >> word)
            {	
              std::istringstream istr(word);

              if (counter%5 == 0)
              {
                if (!(istr >> dummy_int))
                  curv_first.push_back(atoi(vars_[word][0].c_str()));
                else
                  curv_first.push_back(atoi(word.c_str()));
              }
              else if (counter%5 == 1)
              {
                if (!(istr >> dummy_int))
                  curv_second.push_back(atoi(vars_[word][0].c_str()));
                else
                  curv_second.push_back(atoi(word.c_str()));
              }
              else if (counter%5 == 2)
              {
                if (istr >> dummy_dbl)
                {
                  curv_third.push_back(atof(word.c_str()));

                  curv_nurbs.push_back(false);

                  curv_inner_pts.push_back("none");
                  curv_knots.push_back("none");

                  counter += 2;
                }
                else
                {	
                  curv_third.push_back(atof(vars_[word][0].c_str()));
                  curv_nurbs.push_back(true);
                }	
              }
              else if (counter%5 == 3)
                curv_inner_pts.push_back(word);
              else
                curv_knots.push_back(word);

              ++counter;
            }
          }
          else if (isRef)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if (counter%2 == 0)
              {
                if (!(istr >> dummy_int))
                  ref_elt.push_back(atoi(vars_[word][0].c_str()));
                else
                  ref_elt.push_back(atoi(word.c_str()));
              }
              else
              {
                if (!(istr >> dummy_int))
                  ref_type.push_back(atoi(vars_[word][0].c_str()));
                else	
                  ref_type.push_back(atoi(word.c_str()));
              }

              ++counter;
            }
          }
          else if (isVar)
          {
            while (stream >> word)
            {
              vars_[temp_word].push_back(word);
              ++counter;
            }
          }
        }
      }

      assert(x_vertex.size() == y_vertex.size());
      n_vert = x_vertex.size();

      assert(en1.size() == en2.size());
      assert(en2.size() == en3.size());
      assert(en3.size() == en4.size());
      assert(en4.size() == e_mtl.size());
      n_el = en1.size();

      assert(bdy_first.size() == bdy_second.size());
      assert(bdy_first.size() == bdy_type.size());
      n_bdy = bdy_first.size();

      assert(curv_first.size() == curv_second.size());

      n_curv = curv_first.size();

      assert(ref_elt.size() == ref_type.size());
      n_ref = ref_elt.size();
    }
  }
}