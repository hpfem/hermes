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

# include "mesh_data.h"

namespace Hermes
{
  namespace Hermes2D
  {
    MeshData::MeshData(const std::string &mesh_file) : mesh_file_(mesh_file)
    {
    }

    MeshData::~MeshData() {}

    MeshData::MeshData(const MeshData &m) : mesh_file_(m.mesh_file_), n_vert(m.n_vert), n_el(m.n_el), n_bdy(m.n_bdy), n_curv(m.n_curv), n_ref(m.n_ref)
    {
      vars_ = m.vars_;

      x_vertex = m.x_vertex;
      y_vertex = m.y_vertex;

      en1 = m.en1; en2 = m.en2; en3 = m.en3; en4 = m.en4;
      e_mtl = m.e_mtl;

      bdy_first = m.bdy_first; bdy_second = m.bdy_second;
      bdy_type = m.bdy_type;

      curv_first = m.curv_first; curv_second = m.curv_second;
      curv_third = m.curv_third;
      curv_inner_pts = m.curv_inner_pts;
      curv_knots = m.curv_knots;
      curv_nurbs = m.curv_nurbs;

      ref_elt = m.ref_elt;
      ref_type = m.ref_type;
    }

    MeshData& MeshData::operator = (const MeshData &m)
    {
      assert(&m != this);

      mesh_file_ = m.mesh_file_;
      vars_ = m.vars_;

      n_vert = m.n_vert;
      n_el = m.n_el;
      n_bdy = m.n_bdy;
      n_curv = m.n_curv;
      n_ref = m.n_ref;

      x_vertex = m.x_vertex;
      y_vertex = m.y_vertex;

      en1 = m.en1; en2 = m.en2; en3 = m.en3; en4 = m.en4;
      e_mtl = m.e_mtl;

      bdy_first = m.bdy_first; bdy_second = m.bdy_second;
      bdy_type = m.bdy_type;

      curv_first = m.curv_first; curv_second = m.curv_second;
      curv_third = m.curv_third;
      curv_inner_pts = m.curv_inner_pts;
      curv_knots = m.curv_knots;
      curv_nurbs = m.curv_nurbs;

      ref_elt = m.ref_elt;
      ref_type = m.ref_type;

      return *this;
    }

    void MeshData::strip(std::string& str)
    {
      std::string temp;

      // Remove comments
      if(str.find('#') != str.npos)
        str.erase(str.find('#'));

      // Remove brackets, commas and unnecessary tab spaces
      for (size_t i = 0; i < str.length(); i++)
      {
        //if(str[i] != ' ' && str[i] != '\t' && str[i] != '[' && str[i] != ']' && str[i] != '{' && str[i] != '}' && str[i] != '"')
        if(str[i] != '\t' && str[i] != '[' && str[i] != ']' && str[i] != '{' && str[i] != '}')
        {
          if(str[i] == ',' || str[i] == ';')
            temp.append("\t");
          else if(str[i] == '=')
            temp.append("=\t");
          else
            temp.append(1,str[i]);
        }
      }

      str.assign(temp);
      temp.clear();

      // Remove leading whitespaces
      str.erase(0,str.find_first_not_of("\t "));

      // Remove trailing whitespaces
      if(str.find_last_of("\t ") == (str.size() - 1))
        str.erase(str.find_last_not_of("\t ") + 1);

      // Remove unnecessary blank spaces
      for (size_t i = 0; i < str.length(); i++)
      {
        if(str[i] == ' ')
        {
          if(str[i + 1] != ' ' && str[i + 1] != '\t' && str[i + 1] != '=' && str[i - 1] != ' ' && str[i - 1] != '\t')
            temp.append(1,';'); // Meaningful blank spaces are temporarily replaced with ';'
        }
        else
          temp.append(1,str[i]);
      }

      str.assign(temp);
    }

    std::string MeshData::restore(std::string &str)
    {
      std::string temp;

      for (size_t i = 0; i < str.length(); i++)
      {
        if(str[i] != '"')
        {
          if(str[i] == ';')
            temp.append(1,' ');
          else
            temp.append(1,str[i]);
        }
      }

      str.assign(temp);
      return temp;
    }

    void MeshData::parse_mesh(void)
    {
      std::vector<std::string> varlist;

      int dummy_int;
      double dummy_dbl;

      std::ifstream inFile(mesh_file_.c_str());
      std::string line, word, temp_word, next_word;

      int counter(0);
      bool isVert(false), isElt(false), isBdy(false), isCurv(false), isRef(false), isVar(false);

      while (std::getline(inFile,line))
      {
        // Remove all comments, unnecessary blank spaces, commas and paranthesis
        strip(line);

        if(line.find_first_not_of("\t ") != line.npos)
        {
          std::istringstream stream(line);
          stream >> word;

          if(*word.rbegin() == '=')
          {
            word.erase(word.size() - 1);

            if(word == "vertices")
            {
              isVert = true;
              isElt = false; isBdy = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if(word == "elements")
            {
              isElt = true;
              isVert = false; isBdy = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if(word == "boundaries")
            {
              isBdy = true;
              isVert = false; isElt = false; isCurv = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if(word == "curves")
            {
              isCurv = true;
              isVert = false; isElt = false; isBdy = false; isRef = false; isVar = false;
              counter = -1;
            }
            else if(word == "refinements")
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

              temp_word = restore(word);
            }
          }

          if(counter == -1)
            counter = 0;
          else
          {
            if(isVert)
            {
              std::istringstream istr(word);

              if(!(istr >> dummy_dbl))
                x_vertex.push_back(atof(vars_[restore(word)][0].c_str()));
              else
                x_vertex.push_back(atof(word.c_str()));

              ++counter;
            }
            if(isElt)
            {
              std::istringstream istr(word);

              if(!(istr >> dummy_int))
                en1.push_back(atoi(vars_[restore(word)][0].c_str()));
              else
                en1.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if(isBdy)
            {
              std::istringstream istr(word);

              if(!(istr >> dummy_dbl))
                bdy_first.push_back(atof(vars_[restore(word)][0].c_str()));
              else
                bdy_first.push_back(atof(word.c_str()));

              ++counter;
            }
            else if(isCurv)
            {
              std::istringstream istr(word);

              if(!(istr >> dummy_int))
                curv_first.push_back(atoi(vars_[restore(word)][0].c_str()));
              else
                curv_first.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if(isRef)
            {
              std::istringstream istr(word);

              if(!(istr >> dummy_int))
                ref_elt.push_back(atoi(vars_[restore(word)][0].c_str()));
              else
                ref_elt.push_back(atoi(word.c_str()));

              ++counter;
            }
            else if(isVar)
            {
              vars_[temp_word].push_back(restore(word));
              ++counter;
            }
          }

          if(isVert)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if(counter%2 == 0)
              {
                if(!(istr >> dummy_dbl))
                  x_vertex.push_back(atof(vars_[restore(word)][0].c_str()));
                else
                  x_vertex.push_back(atof(word.c_str()));
              }
              else
              {
                if(!(istr >> dummy_dbl))
                  y_vertex.push_back(atof(vars_[restore(word)][0].c_str()));
                else
                  y_vertex.push_back(atof(word.c_str()));
              }

              ++counter;
            }
          }
          else if(isElt)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if(counter%5 == 0)
              {
                if(!(istr >> dummy_int))
                  en1.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  en1.push_back(atoi(word.c_str()));
              }
              else if(counter%5 == 1)
              {
                if(!(istr >> dummy_int))
                  en2.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  en2.push_back(atoi(word.c_str()));
              }
              else if(counter%5 == 2)
              {
                if(!(istr >> dummy_int))
                  en3.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  en3.push_back(atoi(word.c_str()));
              }
              else if(counter%5 == 3)
              {
                if(!(istr >> dummy_int))
                {
                  en4.push_back(-1);
                  e_mtl.push_back(restore(word));

                  ++counter;
                }
                else
                  en4.push_back(atoi(word.c_str()));
              }
              else
                e_mtl.push_back(restore(word));

              ++counter;
            }
          }
          else if(isBdy)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if(counter%3 == 0)
              {
                if(!(istr >> dummy_int))
                  bdy_first.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  bdy_first.push_back(atoi(word.c_str()));
              }
              else if(counter%3 == 1)
              {
                if(!(istr >> dummy_int))
                  bdy_second.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  bdy_second.push_back(atoi(word.c_str()));
              }
              else
                bdy_type.push_back(restore(word));

              ++counter;
            }
          }
          else if(isCurv)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if(counter%5 == 0)
              {
                if(!(istr >> dummy_int))
                  curv_first.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  curv_first.push_back(atoi(word.c_str()));
              }
              else if(counter%5 == 1)
              {
                if(!(istr >> dummy_int))
                  curv_second.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  curv_second.push_back(atoi(word.c_str()));
              }
              else if(counter%5 == 2)
              {
                if(istr >> dummy_dbl)
                {
                  curv_third.push_back(atof(word.c_str()));
                }
                else
                {
                  curv_third.push_back(atof(vars_[restore(word)][0].c_str()));
                }

                stream >> next_word;

                if(next_word == "")
                {
                  curv_nurbs.push_back(false);

                  curv_inner_pts.push_back("none");
                  curv_knots.push_back("none");

                  counter += 2;
                }
                else
                {
                  curv_nurbs.push_back(true);
                }
              }
              else if(counter%5 == 3)
              {
                curv_inner_pts.push_back(restore(next_word));
                curv_knots.push_back(restore(word));
                ++counter;
              }

              ++counter;
            }
          }
          else if(isRef)
          {
            while (stream >> word)
            {
              std::istringstream istr(word);

              if(counter%2 == 0)
              {
                if(!(istr >> dummy_int))
                  ref_elt.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  ref_elt.push_back(atoi(word.c_str()));
              }
              else
              {
                if(!(istr >> dummy_int))
                  ref_type.push_back(atoi(vars_[restore(word)][0].c_str()));
                else
                  ref_type.push_back(atoi(word.c_str()));
              }

              ++counter;
            }
          }
          else if(isVar)
          {
            while (stream >> word)
            {
              vars_[temp_word].push_back(restore(word));
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