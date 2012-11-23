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

#include "filter.h"
#include "forms.h"
#include "quad.h"
#include "refmap.h"
#include "traverse.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Filter<Scalar>::Filter(MeshFunction<Scalar>** solutions, int num) : MeshFunction<Scalar>()
    {
      this->num = num;
      if(num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      for(int i = 0; i < this->num; i++)
        this->sln[i] = solutions[i];
      this->init();
    }

    template<typename Scalar>
    Filter<Scalar>::Filter(const Hermes::vector<MeshFunction<Scalar>*>& solutions) : MeshFunction<Scalar>()
    {
      this->num = solutions.size();
      if(num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      for(int i = 0; i < this->num; i++)
        this->sln[i] = solutions.at(i);
      this->init();
    }

    template<typename Scalar>
    Filter<Scalar>::Filter(const Hermes::vector<Solution<Scalar>*>& solutions) : MeshFunction<Scalar>()
    {
      this->num = solutions.size();
      if(num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      for(int i = 0; i < this->num; i++)
        this->sln[i] = solutions.at(i);
      this->init();
    }

    template<typename Scalar>
    void Filter<Scalar>::init(const Hermes::vector<MeshFunction<Scalar>*>& solutions)
    {
      this->num = solutions.size();
      if(num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      for(int i = 0; i < this->num; i++)
        this->sln[i] = solutions.at(i);
      this->init();
    }

    template<typename Scalar>
    void Filter<Scalar>::setDeleteSolutions()
    {
      this->deleteSolutions = true;
    }

    template<typename Scalar>
    void Filter<Scalar>::init()
    {
      // construct the union mesh, if necessary
      const Mesh* meshes[10];
      for(int i = 0; i < this->num; i++)
        meshes[i] = this->sln[i]->get_mesh();
      this->mesh = meshes[0];
      unimesh = false;

      for (int i = 1; i < num; i++)

      {
        if(meshes[i] == NULL)
        {
          this->warn("You may be initializing a Filter with Solution that is missing a Mesh.");
          throw Hermes::Exceptions::Exception("this->meshes[%d] is NULL in Filter<Scalar>::init().", i);
        }
        if(meshes[i]->get_seq() != this->mesh->get_seq())
        {
          unimesh = true;
          break;
        }
      }

      if(unimesh)
      {
        Traverse trav;
        trav.begin(num, meshes);
        this->mesh = new Mesh;
        unidata = trav.construct_union_mesh(const_cast<Hermes2D::Mesh*>(this->mesh));
        trav.finish();
      }

      // misc init
      this->num_components = 1;
      this->order = 0;

      memset(sln_sub, 0, sizeof(sln_sub));
      set_quad_2d(&g_quad_2d_std);

      this->deleteSolutions = false;
    }

    template<typename Scalar>
    Filter<Scalar>::~Filter()
    {
      free();
      if(this->deleteSolutions)
      {
        for(int i = 0; i < this->num; i++)
          delete this->sln[i];
      }
    }

    template<typename Scalar>
    void Filter<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      MeshFunction<Scalar>::set_quad_2d(quad_2d);
      for (int i = 0; i < num; i++)
        this->sln[i]->set_quad_2d(quad_2d); // nodup
    }

    template<typename Scalar>
    void Filter<Scalar>::set_active_element(Element* e)
    {
      MeshFunction<Scalar>::set_active_element(e);
      if(!unimesh)
      {
        for (int i = 0; i < num; i++)
          this->sln[i]->set_active_element(e); // nodup
        memset(sln_sub, 0, sizeof(sln_sub));
      }
      else
      {
        for (int i = 0; i < num; i++)
        {
          this->sln[i]->set_active_element(unidata[i][e->id].e);
          this->sln[i]->set_transform(unidata[i][e->id].idx);
          sln_sub[i] = this->sln[i]->get_transform();
        }
      }

      for(typename std::map<uint64_t, LightArray<struct Filter<Scalar>::Node*>*>::iterator it = tables[this->cur_quad].begin(); it != tables[this->cur_quad].end(); it++)
      {
        for(unsigned int l = 0; l < it->second->get_size(); l++)
          if(it->second->present(l))
            ::free(it->second->get(l));
        delete it->second;
      }
      tables[this->cur_quad].clear();

      this->sub_tables = &tables[this->cur_quad];
      this->update_nodes_ptr();

      this->order = 20; // fixme
    }

    template<typename Scalar>
    void Filter<Scalar>::free()
    {
      for (int i = 0; i < num; i++)
      {
        for(typename std::map<uint64_t, LightArray<struct Filter<Scalar>::Node*>*>::iterator it = tables[i].begin(); it != tables[i].end(); it++)
        {
          for(unsigned int l = 0; l < it->second->get_size(); l++)
            if(it->second->present(l))
              ::free(it->second->get(l));
          delete it->second;
        }
        tables[i].clear();
      }

      if(unimesh)
      {
        delete this->mesh;
        for (int i = 0; i < num; i++)
          ::free(unidata[i]);
        delete [] unidata;
      }
    }

    template<typename Scalar>
    void Filter<Scalar>::reinit()
    {
      free();
      init();
    }

    template<typename Scalar>
    void Filter<Scalar>::push_transform(int son)
    {
      MeshFunction<Scalar>::push_transform(son);
      for (int i = 0; i < num; i++)
      {
        // sln_sub[i] contains the value sln[i]->sub_idx, which the Filter thinks
        // the solution has, or at least had last time. If the filter graph is
        // cyclic, it could happen that some solutions would get pushed the transform
        // more than once. This mechanism prevents it. If the filter sees that the
        // solution already has a different sub_idx than it thinks it should have,
        // it assumes someone else has already pushed the correct transform. This
        // is actually the case for cyclic filter graphs and filters used in multi-mesh
        // computation.

        if(this->sln[i]->get_transform() == sln_sub[i])
          this->sln[i]->push_transform(son);
        sln_sub[i] = this->sln[i]->get_transform();
      }
    }

    template<typename Scalar>
    void Filter<Scalar>::pop_transform()
    {
      MeshFunction<Scalar>::pop_transform();
      for (int i = 0; i < num; i++)
      {
        if(this->sln[i]->get_transform() == sln_sub[i] && sln_sub[i] != 0)
          this->sln[i]->pop_transform();
        sln_sub[i] = this->sln[i]->get_transform();
      }
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::SimpleFilter(const Hermes::vector<MeshFunction<Scalar>*>& solutions, const Hermes::vector<int>& items)
    {
      this->num = solutions.size();
      if(this->num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      if(items.size() != (unsigned) this->num)
        if(items.size() > 0)
          throw Hermes::Exceptions::Exception("Attempt to create an instance of SimpleFilter with different supplied number of MeshFunctions than the number of types of data used from them.");

      for(int i = 0; i < this->num; i++)

      {
        this->sln[i] = solutions.at(i);
        if(items.size() > 0)
          this->item[i] = items.at(i);
        else
          this->item[i] =  H2D_FN_VAL;
      }
      this->init();
      init_components();
      this->deleteSolutions = false;
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::SimpleFilter(const Hermes::vector<Solution<Scalar>*>& solutions, const Hermes::vector<int>& items)
    {
      this->num = solutions.size();
      if(this->num > 10)
        throw Hermes::Exceptions::Exception("Attempt to create an instance of Filter with more than 10 MeshFunctions.");
      if(items.size() != (unsigned) this->num)
        if(items.size() > 0)
          throw Hermes::Exceptions::Exception("Attempt to create an instance of SimpleFilter with different supplied number of MeshFunctions than the number of types of data used from them.");

      for(int i = 0; i < this->num; i++)
      {
        this->sln[i] = solutions.at(i);
        if(items.size() > 0)
          this->item[i] = items.at(i);
        else
          this->item[i] = H2D_FN_VAL;
      }
      this->init();
      init_components();
      this->deleteSolutions = false;
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::~SimpleFilter()
    {
    }

    template<typename Scalar>
    void SimpleFilter<Scalar>::init_components()
    {
      bool vec1 = false, vec2 = false;
      for (int i = 0; i < this->num; i++)
      {
        if(this->sln[i]->get_num_components() > 1) vec1 = true;
        if((item[i] & H2D_FN_COMPONENT_0) && (item[i] & H2D_FN_COMPONENT_1)) vec2 = true;
        if(this->sln[i]->get_num_components() == 1) item[i] &= H2D_FN_COMPONENT_0;
      }
      this->num_components = (vec1 && vec2) ? 2 : 1;
    }

    template<typename Scalar>
    void SimpleFilter<Scalar>::precalculate(int order, int mask)
    {
      if(mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
        throw Hermes::Exceptions::Exception("Filter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->element->get_mode());
      struct Function<Scalar>::Node* node = this->new_node(H2D_FN_VAL, np);

      // precalculate all solutions
      for (int i = 0; i < this->num; i++)
        this->sln[i]->set_quad_order(order, item[i]);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain corresponding tables
        Scalar* tab[10];
        for (int i = 0; i < this->num; i++)
        {
          int a = 0, b = 0, mask = item[i];
          if(mask >= 0x40) { a = 1; mask >>= 6; }
          while (!(mask & 1)) { mask >>= 1; b++; }
          tab[i] = this->sln[i]->get_values(this->num_components == 1 ? a : j, b);
          if(tab[i] == NULL) throw Hermes::Exceptions::Exception("Value of 'item%d' is incorrect in filter definition.", i + 1);
        }

        Hermes::vector<Scalar*> values;
        for(int i = 0; i < this->num; i++)
          values.push_back(tab[i]);

        // apply the filter
        filter_fn(np, values, node->values[j][0]);
      }

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      this->cur_node = node;
    }

    template<typename Scalar>
    Func<Scalar>* SimpleFilter<Scalar>::get_pt_value(double x, double y)
    {
      Scalar val[10];
      for (int i = 0; i < this->num; i++)
        val[i] = this->sln[i]->get_pt_value(x, y)->val[0];

      Func<Scalar>* toReturn = new Func<Scalar>(1, 1);

      Scalar result;
      Hermes::vector<Scalar*> values;
      for(int i = 0; i < this->num; i++)
        values.push_back(&val[i]);

      // apply the filter
      filter_fn(1, values, &result);

      toReturn->val[0] = result;
      return toReturn;
    }

    template<typename Scalar>
    DXFilter<Scalar>::DXFilter(const Hermes::vector<MeshFunction<Scalar>*>& solutions) : DXDYFilter<Scalar>(solutions)
    {
    }

    template<typename Scalar>
    void DXFilter<Scalar>::filter_fn(int n, Hermes::vector<Scalar *> values, Hermes::vector<Scalar *> dx, Hermes::vector<Scalar *> dy, Scalar* rslt, Scalar* rslt_dx, Scalar* rslt_dy)
    {
      for (int i = 0; i < n; i++)
      {
        rslt_dx[i] = 0;
        rslt_dy[i] = 0;
        for(unsigned int j = 0; j < values.size(); j++)
        {
          rslt_dx[i] += dx.at(j)[i];
          rslt_dy[i] += dy.at(j)[i];
        }
      }
    }

    template<typename Scalar>
    MeshFunction<Scalar>* DXFilter<Scalar>::clone() const
    {
      Hermes::vector<MeshFunction<Scalar>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
      }
      DXFilter<Scalar>* filter = new DXFilter<Scalar>(slns);
      filter->setDeleteSolutions();
      return filter;
    }


    ComplexFilter::ComplexFilter(MeshFunction<std::complex<double> >* solution, int item) : Filter<double>()
    {
      this->num = 0;
      this->unimesh = false;
      this->sln_complex = solution;
      this->num_components = solution->get_num_components();
      this->mesh = solution->get_mesh();
      set_quad_2d(&g_quad_2d_std);
    }

    void ComplexFilter::free()
    {
      for(std::map<uint64_t, LightArray<struct Filter<double>::Node*>*>::iterator it = tables[this->cur_quad].begin(); it != tables[this->cur_quad].end(); it++)
      {
        for(unsigned int l = 0; l < it->second->get_size(); l++)
          if(it->second->present(l))
            ::free(it->second->get(l));
        delete it->second;
      }
      tables[this->cur_quad].clear();
      if(this->deleteSolutions)
        delete this->sln_complex;
    }

    void ComplexFilter::set_quad_2d(Quad2D* quad_2d)
    {
      MeshFunction<double>::set_quad_2d(quad_2d);
      this->sln_complex->set_quad_2d(quad_2d);
    }

    void ComplexFilter::set_active_element(Element* e)
    {
      MeshFunction<double>::set_active_element(e);

      this->sln_complex->set_active_element(e);

      memset(sln_sub, 0, sizeof(sln_sub));

      for(std::map<uint64_t, LightArray<struct Filter<double>::Node*>*>::iterator it = tables[this->cur_quad].begin(); it != tables[this->cur_quad].end(); it++)
      {
        for(unsigned int l = 0; l < it->second->get_size(); l++)
          if(it->second->present(l))
            ::free(it->second->get(l));
        delete it->second;
      }
      tables[this->cur_quad].clear();

      this->sub_tables = &tables[this->cur_quad];

      this->update_nodes_ptr();

      this->order = 20; // fixme
    }

    void ComplexFilter::push_transform(int son)
    {
      MeshFunction<double>::push_transform(son);
      this->sln_complex->push_transform(son);
    }

    void ComplexFilter::pop_transform()
    {
      MeshFunction<double>::pop_transform();
      this->sln_complex->pop_transform();
    }

    void ComplexFilter::precalculate(int order, int mask)
    {
      if(mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
        throw Hermes::Exceptions::Exception("Filter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->element->get_mode());
      struct Function<double>::Node* node = this->new_node(H2D_FN_VAL, np);

      this->sln_complex->set_quad_order(order, H2D_FN_VAL);

      // obtain corresponding tables
      filter_fn(np, this->sln_complex->get_values(0, 0), node->values[0][0]);
      if(num_components > 1)
        filter_fn(np, this->sln_complex->get_values(1, 0), node->values[1][0]);

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }

      this->cur_node = node;
    }

    Func<double>* ComplexFilter::get_pt_value(double x, double y)
    {
      Func<std::complex<double> >* val = this->sln_complex->get_pt_value(x, y);

      Func<double>* toReturn = new Func<double>(1, this->num_components);

      double result;

      // apply the filter
      filter_fn(1, val->val, &result);

      if(this->num_components == 1)
      {
        toReturn->val[0] = result;
        toReturn->dx[0] = 0.0;
        toReturn->dy[0] = 0.0;
      }
      else
      {
        this->warn("ComplexFilter only implemented for scalar functions."); 
      }
      return toReturn;
    }

    template<typename Scalar>
    DXDYFilter<Scalar>::DXDYFilter(const Hermes::vector<MeshFunction<Scalar>*>& solutions) : Filter<Scalar>(solutions)
    {
      init_components();
    }

    template<typename Scalar>
    DXDYFilter<Scalar>::DXDYFilter(const Hermes::vector<Solution<Scalar>*>& solutions) : Filter<Scalar>(solutions)
    {
      init_components();
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::init_components()
    {
      this->num_components = this->sln[0]->get_num_components();
      for (int i = 1; i < this->num; i++)
        if(this->sln[i]->get_num_components() != this->num_components)
          throw Hermes::Exceptions::Exception("Filter: Solutions do not have the same number of components!");
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::init(const Hermes::vector<MeshFunction<Scalar>*>& solutions)
    {
      Filter<Scalar>::init(solutions);
      init_components();
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::precalculate(int order, int mask)
    {
      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->element->get_mode());
      struct Function<Scalar>::Node* node = this->new_node(H2D_FN_DEFAULT, np);

      // precalculate all solutions
      for (int i = 0; i < this->num; i++)
        this->sln[i]->set_quad_order(order, H2D_FN_DEFAULT);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain solution tables
        Scalar *val[10], *dx[10], *dy[10];
        for (int i = 0; i < this->num; i++)
        {
          val[i] = this->sln[i]->get_fn_values(j);
          dx[i]  = this->sln[i]->get_dx_values(j);
          dy[i]  = this->sln[i]->get_dy_values(j);
        }

        Hermes::vector<Scalar *> values_vector;
        Hermes::vector<Scalar *> dx_vector;
        Hermes::vector<Scalar *> dy_vector;

        for(int i = 0; i < this->num; i++)

        {
          values_vector.push_back(val[i]);
          dx_vector.push_back(dx[i]);
          dy_vector.push_back(dy[i]);
        }

        // apply the filter
        filter_fn(np, values_vector, dx_vector, dy_vector, node->values[j][0], node->values[j][1], node->values[j][2]);
      }

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      this->cur_node = node;
    }

    template<typename Scalar>
    Func<Scalar>* DXDYFilter<Scalar>::get_pt_value(double x, double y)
    {
      return 0;
    }

    template<typename Scalar>
    void MagFilter<Scalar>::filter_fn(int n, Hermes::vector<Scalar*> values, Scalar* result)
    {
      for (int i = 0; i < n; i++)

      {
        result[i] = 0;
        for(unsigned int j = 0; j < values.size(); j++)
          result[i] += sqr(values.at(j)[i]);
        result[i] = sqrt(result[i]);
      }
    };

    template<typename Scalar>
    MagFilter<Scalar>::MagFilter(Hermes::vector<MeshFunction<Scalar>*> solutions, Hermes::vector<int> items) : SimpleFilter<Scalar>(solutions, items)
    {
    };

    template<typename Scalar>
    MagFilter<Scalar>::MagFilter(MeshFunction<Scalar>* sln1, int item1)
      : SimpleFilter<Scalar>(Hermes::vector<MeshFunction<Scalar>*>(sln1, sln1),
      Hermes::vector<int>(item1 & H2D_FN_COMPONENT_0, item1 & H2D_FN_COMPONENT_1))
    {
      if(sln1->get_num_components() < 2)
        throw Hermes::Exceptions::Exception("The single-argument constructor is intended for vector-valued solutions.");
    };

    template<typename Scalar>
    MeshFunction<Scalar>* MagFilter<Scalar>::clone() const
    {
      Hermes::vector<MeshFunction<Scalar>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      MagFilter<Scalar>* filter = new MagFilter<Scalar>(slns, items);
      filter->setDeleteSolutions();
      return filter;
    }

    void TopValFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for(unsigned int j = 0; j < values.size(); j++)
          if(values.at(j)[i] > limits[j])
            result[i] = limits[j];
          else
            result[i] = values.at(j)[i];
      }
    };

    TopValFilter::TopValFilter(Hermes::vector<MeshFunction<double>*> solutions, Hermes::vector<double> limits, Hermes::vector<int> items) : SimpleFilter<double>(solutions, items), limits(limits)
    {
    };

    TopValFilter::TopValFilter(MeshFunction<double>* sln, double limit, int item)
      : SimpleFilter<double>()
    {
      this->limits.push_back(limit);
      this->sln[0] = sln;
      this->item[0] = item;
      this->num = 1;
      Filter<double>::init();
    };

    MeshFunction<double>* TopValFilter::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      TopValFilter* filter = new TopValFilter(slns, limits, items);
      filter->setDeleteSolutions();
      return filter;
    }

    void BottomValFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for(unsigned int j = 0; j < values.size(); j++)
          if(values.at(j)[i] < limits[j])
            result[i] = limits[j];
          else
            result[i] = values.at(j)[i];
      }
    };

    BottomValFilter::BottomValFilter(Hermes::vector<MeshFunction<double>*> solutions, Hermes::vector<double> limits, Hermes::vector<int> items) : SimpleFilter<double>(solutions, items), limits(limits)
    {
    };

    BottomValFilter::BottomValFilter(MeshFunction<double>* sln, double limit, int item)
      : SimpleFilter<double>()
    {
      this->limits.push_back(limit);
      this->sln[0] = sln;
      this->item[0] = item;
      this->num = 1;
      Filter<double>::init();
    };

    MeshFunction<double>* BottomValFilter::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      BottomValFilter* filter = new BottomValFilter(slns, limits, items);
      filter->setDeleteSolutions();
      return filter;
    }

    void ValFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for(unsigned int j = 0; j < values.size(); j++)
          if(values.at(j)[i] < low_limits[j])
            result[i] = low_limits[j];
          else
            if(values.at(j)[i] > high_limits[j])
              result[i] = high_limits[j];
            else
              result[i] = values.at(j)[i];
      }
    };

    ValFilter::ValFilter(Hermes::vector<MeshFunction<double>*> solutions, Hermes::vector<double> low_limits, Hermes::vector<double> high_limits, Hermes::vector<int> items) : SimpleFilter<double>(solutions, items), low_limits(low_limits), high_limits(high_limits)
    {
    };

    ValFilter::ValFilter(MeshFunction<double>* sln, double low_limit, double high_limit, int item)
      : SimpleFilter<double>()
    {
      this->low_limits.push_back(low_limit);
      this->high_limits.push_back(high_limit);
      this->sln[0] = sln;
      this->item[0] = item;
      this->num = 1;
      Filter<double>::init();
    };

    MeshFunction<double>* ValFilter::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      ValFilter* filter = new ValFilter(slns, low_limits, high_limits, items);
      filter->setDeleteSolutions();
      return filter;
    }

    template<typename Scalar>
    void DiffFilter<Scalar>::filter_fn(int n, Hermes::vector<Scalar*> values, Scalar* result)
    {
      for (int i = 0; i < n; i++) result[i] = values.at(0)[i] - values.at(1)[i];
    };

    template<typename Scalar>
    DiffFilter<Scalar>::DiffFilter(Hermes::vector<MeshFunction<Scalar>*> solutions, Hermes::vector<int> items) : SimpleFilter<Scalar>(solutions, items) {}

    template<typename Scalar>
    MeshFunction<Scalar>* DiffFilter<Scalar>::clone() const
    {
      Hermes::vector<MeshFunction<Scalar>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      DiffFilter* filter = new DiffFilter<Scalar>(slns, items);
      filter->setDeleteSolutions();
      return filter;
    }

    template<typename Scalar>
    void SumFilter<Scalar>::filter_fn(int n, Hermes::vector<Scalar*> values, Scalar* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          result[i] += values.at(j)[i];
      }
    };

    template<typename Scalar>
    SumFilter<Scalar>::SumFilter(Hermes::vector<MeshFunction<Scalar>*> solutions, Hermes::vector<int> items) : SimpleFilter<Scalar>(solutions, items) {}

    template<typename Scalar>
    MeshFunction<Scalar>* SumFilter<Scalar>::clone() const
    {
      Hermes::vector<MeshFunction<Scalar>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      SumFilter<Scalar>* filter = new SumFilter<Scalar>(slns, items);
      filter->setDeleteSolutions();
      return filter;
    }

    template<>
    void SquareFilter<double>::filter_fn(int n, Hermes::vector<double *> v1, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = sqr(v1.at(0)[i]);
    };

    template<>
    void SquareFilter<std::complex<double> >::filter_fn(int n, Hermes::vector<std::complex<double> *> v1, std::complex<double> * result)
    {
      for (int i = 0; i < n; i++)
        result[i] = std::norm(v1.at(0)[i]);
    };

    template<typename Scalar>
    SquareFilter<Scalar>::SquareFilter(Hermes::vector<MeshFunction<Scalar>*> solutions, Hermes::vector<int> items)
      : SimpleFilter<Scalar>(solutions, items)
    {
      if(solutions.size() > 1)
        throw Hermes::Exceptions::Exception("SquareFilter only supports one MeshFunction.");
    };

    template<typename Scalar>
    MeshFunction<Scalar>* SquareFilter<Scalar>::clone() const
    {
      Hermes::vector<MeshFunction<Scalar>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      SquareFilter<Scalar>* filter = new SquareFilter<Scalar>(slns, items);
      filter->setDeleteSolutions();
      return filter;
    }

    void AbsFilter::filter_fn(int n, Hermes::vector<double*> v1, double * result)
    {
      for (int i = 0; i < n; i++)
        result[i] = std::abs(v1.at(0)[i]);
    };

    AbsFilter::AbsFilter(Hermes::vector<MeshFunction<double>*> solutions, Hermes::vector<int> items)
      : SimpleFilter<double>(solutions, items)
    {
      if(solutions.size() > 1)
        throw Hermes::Exceptions::Exception("AbsFilter only supports one MeshFunction.");
    };

    AbsFilter::AbsFilter(MeshFunction<double>* solution)
      : SimpleFilter<double>()
    {
      this->num = 1;
      this->sln[0] = solution;

      this->item[0] =  H2D_FN_VAL;

      this->init();
      init_components();
    };

    MeshFunction<double>* AbsFilter::clone() const
    {
      Hermes::vector<MeshFunction<double>*> slns;
      Hermes::vector<int> items;
      for(int i = 0; i < this->num; i++)
      {
        slns.push_back(this->sln[i]->clone());
        items.push_back(this->item[i]);
      }
      AbsFilter* filter = new AbsFilter(slns, items);
      filter->setDeleteSolutions();
      return filter;
    }

    void RealFilter::filter_fn(int n, std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = values[i].real();
    };

    MeshFunction<double>* RealFilter::clone() const
    {
      RealFilter* filter = new RealFilter(this->sln_complex->clone(), this->item);
      filter->setDeleteSolutions();
      return filter;
    }

    RealFilter::RealFilter(MeshFunction<std::complex<double> >* solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    void ImagFilter::filter_fn(int n, std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = values[i].imag();
    };

    ImagFilter::ImagFilter(MeshFunction<std::complex<double> >* solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    MeshFunction<double>* ImagFilter::clone() const
    {
      ImagFilter* filter = new ImagFilter(this->sln_complex->clone(), this->item);
      filter->setDeleteSolutions();
      return filter;
    }

    void ComplexAbsFilter::filter_fn(int n, std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = sqrt(sqr(values[i].real()) + sqr(values[i].imag()));
    };

    MeshFunction<double>* ComplexAbsFilter::clone() const
    {
      ComplexAbsFilter* filter = new ComplexAbsFilter(this->sln_complex->clone(), this->item);
      filter->setDeleteSolutions();
      return filter;
    }

    ComplexAbsFilter::ComplexAbsFilter(MeshFunction<std::complex<double> >* solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    void AngleFilter::filter_fn(int n, Hermes::vector<std::complex<double>*> v1, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = atan2( v1.at(0)[i].imag(), v1.at(0)[i].real() );
    };

    AngleFilter::AngleFilter(Hermes::vector<MeshFunction<std::complex<double> >*> solutions, Hermes::vector<int> items)
      : SimpleFilter<std::complex<double> >(solutions, items)
    {
      if(solutions.size() > 1)
        throw Hermes::Exceptions::Exception("RealFilter only supports one MeshFunction.");
    };

    void VonMisesFilter::precalculate(int order, int mask)
    {
      if(mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
        throw Hermes::Exceptions::Exception("VonMisesFilter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->element->get_mode());
      Filter<double>::Node* node = new_node(H2D_FN_VAL_0, np);

      this->sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
      this->sln[1]->set_quad_order(order, H2D_FN_DX | H2D_FN_DY);

      double *dudx, *dudy, *dvdx, *dvdy;
      this->sln[0]->get_dx_dy_values(dudx, dudy);
      this->sln[1]->get_dx_dy_values(dvdx, dvdy);
      double *uval = this->sln[0]->get_fn_values();
      update_refmap();
      double *x = refmap->get_phys_x(order);

      for (int i = 0; i < np; i++)
      {
        // stress tensor
        double tz = lambda*(dudx[i] + dvdy[i]);
        double tx = tz + 2*mu*dudx[i];
        double ty = tz + 2*mu*dvdy[i];
        if(cyl) tz += 2*mu*uval[i] / x[i];
        double txy = mu*(dudy[i] + dvdx[i]);

        // Von Mises stress
        node->values[0][0][i] = 1.0/sqrt(2.0) * sqrt(sqr(tx - ty) + sqr(ty - tz) + sqr(tz - tx) + 6*sqr(txy));
      }

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      cur_node = node;
    }

    Func<double>* VonMisesFilter::get_pt_value(double x, double y)
    {
      return 0;
    }

    VonMisesFilter::VonMisesFilter(Hermes::vector<MeshFunction<double>*> solutions, double lambda, double mu,
      int cyl, int item1, int item2)
      : Filter<double>(solutions)
    {
      this->mu = mu;
      this->lambda = lambda;
      this->cyl = cyl;
      this->item1 = item1;
      this->item2 = item2;
    }

    VonMisesFilter::VonMisesFilter(MeshFunction<double>** solutions, int num, double lambda, double mu,
        int cyl, int item1, int item2): Filter<double>(solutions, num)
    {
      this->mu = mu;
      this->lambda = lambda;
      this->cyl = cyl;
      this->item1 = item1;
      this->item2 = item2;
    }

    MeshFunction<double>* VonMisesFilter::clone() const
    {
      MeshFunction<double>** slns = new MeshFunction<double>*[num];
      for(int i = 0; i < num; i++)
        slns[i] = sln[i]->clone();
      VonMisesFilter* filter = new VonMisesFilter(slns, num, lambda, mu, cyl, item1, item2);
      filter->setDeleteSolutions();
      return filter;
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::precalculate(int order, int mask)
    {
      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->element->get_mode());
      struct Filter<Scalar>::Node* node = this->new_node(H2D_FN_DEFAULT, np);

      // precalculate all solutions
      for (int i = 0; i < this->num; i++)
        this->sln[i]->set_quad_order(order);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain solution tables
        Scalar *val[4], *dx[4], *dy[4];
        for (int i = 0; i < this->num; i++)
        {
          val[i] = this->sln[i]->get_fn_values(j);
          dx[i]  = this->sln[i]->get_dx_values(j);
          dy[i]  = this->sln[i]->get_dy_values(j);
        }
        if(this->num == 2)
          for (int i = 0; i < np; i++)
          {
            node->values[j][0][i] = tau_frac * (val[1][i] - val[0][i]) + val[1][i];
            node->values[j][1][i] = tau_frac * (dx[1][i]  - dx[0][i])  + dx[1][i];
            node->values[j][2][i] = tau_frac * (dy[1][i]  - dy[0][i])  + dy[1][i];
          }
        else
          for (int i = 0; i < np; i++)
          {
            node->values[j][0][i] = val[0][i];
            node->values[j][1][i] = dx[0][i];
            node->values[j][2][i] = dy[0][i];
          }
      }

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      this->cur_node = node;
    }

    template<typename Scalar>
    Func<Scalar>* LinearFilter<Scalar>::get_pt_value(double x, double y)
    {
      return 0;
    }

    template<typename Scalar>
    LinearFilter<Scalar>::LinearFilter(MeshFunction<Scalar>* old) // : Filter(old)
    {
      init_components();
    }

    template<typename Scalar>
    LinearFilter<Scalar>::LinearFilter(MeshFunction<Scalar>* older, MeshFunction<Scalar>* old, double tau_frac)
      : Filter<Scalar>(Hermes::vector<MeshFunction<Scalar>*>(older, old))
    {
      this->tau_frac = tau_frac;
      init_components();
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::init_components()
    {
      this->num_components = this->sln[0]->get_num_components();
      for (int i = 1; i < this->num; i++)
        if(this->sln[i]->get_num_components() != this->num_components)
          throw Hermes::Exceptions::Exception("Filter: Solutions do not have the same number of components!");
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::set_active_element(Element* e)
    {
      Filter<Scalar>::set_active_element(e);

      this->order = 0;
      for (int i = 0; i < this->num; i++)
      {
        int o = this->sln[i]->get_fn_order();
        if(o > this->order) this->order = o;
      }
    }

    template class HERMES_API Filter<double>;
    template class HERMES_API Filter<std::complex<double> >;
    template class HERMES_API SimpleFilter<double>;
    template class HERMES_API SimpleFilter<std::complex<double> >;
    template class HERMES_API DXDYFilter<double>;
    template class HERMES_API DXDYFilter<std::complex<double> >;
    template class HERMES_API DXFilter<double>;
    template class HERMES_API DXFilter<std::complex<double> >;
    template class HERMES_API MagFilter<double>;
    template class HERMES_API MagFilter<std::complex<double> >;
    template class HERMES_API DiffFilter<double>;
    template class HERMES_API DiffFilter<std::complex<double> >;
    template class HERMES_API SumFilter<double>;
    template class HERMES_API SumFilter<std::complex<double> >;
    template class HERMES_API SquareFilter<double>;
    template class HERMES_API SquareFilter<std::complex<double> >;
  }
}