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
    Filter<Scalar>::Filter()
    {
    }

    template<typename Scalar>
    Filter<Scalar>::Filter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions) : MeshFunction<Scalar>(), solutions(solutions)
    {
      this->init();
    }

    template<typename Scalar>
    void Filter<Scalar>::init()
    {
      // construct the union mesh, if necessary
      std::vector<MeshSharedPtr> meshes;
      for (int i = 0; i < this->solutions.size(); i++)
        meshes.push_back(this->solutions[i]->get_mesh());
      this->mesh = meshes[0];

      for (int i = 0; i < this->solutions.size(); i++)
        this->solutions_sub_idx.push_back(0);

      Solution<Scalar>* sln = dynamic_cast<Solution<Scalar>*>(this->solutions[0].get());
      if (sln == nullptr)
        this->space_type = HERMES_INVALID_SPACE;
      else
        this->space_type = sln->get_space_type();

      unimesh = false;

      for (int i = 1; i < this->solutions.size(); i++)
      {
        if (meshes[i] == nullptr)
        {
          this->warn("You may be initializing a Filter with Solution that is missing a Mesh.");
          throw Hermes::Exceptions::Exception("this->meshes[%d] is nullptr in Filter<Scalar>::init().", i);
        }
        if (meshes[i]->get_seq() != this->mesh->get_seq())
        {
          unimesh = true;
          break;
        }

        sln = dynamic_cast<Solution<Scalar>*>(this->solutions[i].get());
        if (sln == nullptr || sln->get_space_type() != this->space_type)
          this->space_type = HERMES_INVALID_SPACE;
      }

      if (unimesh)
      {
        this->mesh = MeshSharedPtr(new Mesh);
        this->unidata = Traverse::construct_union_mesh(this->solutions.size(), &meshes[0], this->mesh);
      }

      // misc init
      this->num_components = 1;
      this->order = 0;

      set_quad_2d(&g_quad_2d_std);
    }

    template<typename Scalar>
    Filter<Scalar>::~Filter()
    {
      this->free();
    }

    template<typename Scalar>
    void Filter<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      MeshFunction<Scalar>::set_quad_2d(quad_2d);
      for (int i = 0; i < this->solutions.size(); i++)
        // nodup
        this->solutions[i]->set_quad_2d(quad_2d);
    }

    template<typename Scalar>
    void Filter<Scalar>::set_active_element(Element* e)
    {
      MeshFunction<Scalar>::set_active_element(e);
      if (!unimesh)
      {
        for (int i = 0; i < this->solutions.size(); i++)
          // nodup
          this->solutions[i]->set_active_element(e);

        for (int i = 0; i < this->solutions.size(); i++)
          solutions_sub_idx[i] = 0;
      }
      else
      {
        for (int i = 0; i < this->solutions.size(); i++)
        {
          this->solutions[i]->set_active_element(unidata[i][e->id].e);
          this->solutions[i]->set_transform(unidata[i][e->id].idx);
          solutions_sub_idx[i] = this->solutions[i]->get_transform();
        }
      }

      // fixme
      this->order = 20;
    }

    template<typename Scalar>
    void Filter<Scalar>::free()
    {
      if (unimesh)
      {
        for (int i = 0; i < this->solutions.size(); i++)
          free_with_check(unidata[i]);
        free_with_check(unidata);
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
      for (int i = 0; i < this->solutions.size(); i++)
      {
        // solutions_sub_idx[i] contains the value sln[i]->sub_idx, which the Filter thinks
        // the solution has, or at least had last time. If the filter graph is
        // cyclic, it could happen that some solutions would get pushed the transform
        // more than once. This mechanism prevents it. If the filter sees that the
        // solution already has a different sub_idx than it thinks it should have,
        // it assumes someone else has already pushed the correct transform. This
        // is actually the case for cyclic filter graphs and filters used in multi-mesh
        // computation.

        if (this->solutions[i]->get_transform() == solutions_sub_idx[i])
          this->solutions[i]->push_transform(son);
        solutions_sub_idx[i] = this->solutions[i]->get_transform();
      }
    }

    template<typename Scalar>
    void Filter<Scalar>::pop_transform()
    {
      MeshFunction<Scalar>::pop_transform();
      for (int i = 0; i < this->solutions.size(); i++)
      {
        if (this->solutions[i]->get_transform() == solutions_sub_idx[i] && solutions_sub_idx[i] != 0)
          this->solutions[i]->pop_transform();
        solutions_sub_idx[i] = this->solutions[i]->get_transform();
      }
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::SimpleFilter() : Filter<Scalar>()
    {
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::SimpleFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions, const std::vector<int> items) : Filter<Scalar>(solutions), items(items)
    {
      if (this->items.size() > 0)
        Hermes::Helpers::check_length(solutions, items);
      else
        for (int i = 0; i < this->solutions.size(); i++)
          this->items.push_back(H2D_FN_VAL);

      this->init();
      this->init_components();
    }

    template<typename Scalar>
    SimpleFilter<Scalar>::~SimpleFilter()
    {
    }

    template<typename Scalar>
    void SimpleFilter<Scalar>::init_components()
    {
      bool vec1 = false, vec2 = false;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        if (this->solutions[i]->get_num_components() > 1)
          vec1 = true;
        
        if ((this->items[i] & H2D_FN_COMPONENT_0) && (this->items[i] & H2D_FN_COMPONENT_1))
          vec2 = true;
        
        if (this->solutions[i]->get_num_components() == 1)
          this->items[i] &= H2D_FN_COMPONENT_0;
      }
      this->num_components = (vec1 && vec2) ? 2 : 1;
    }

    template<typename Scalar>
    void SimpleFilter<Scalar>::precalculate(unsigned short order, unsigned short mask)
    {
#ifdef H2D_USE_SECOND_DERIVATIVES
      if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
#else
      if (mask & (H2D_FN_DX | H2D_FN_DY))
#endif
        throw Hermes::Exceptions::Exception("SimpleFilter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      unsigned char np = quad->get_num_points(order, this->element->get_mode());

      // precalculate all solutions
      for (int i = 0; i < this->solutions.size(); i++)
        this->solutions[i]->set_quad_order(order, this->items[i]);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain corresponding tables
        std::vector<const Scalar*> tab;
        for (int i = 0; i < this->solutions.size(); i++)
        {
          int a = 0, b = 0, mask = this->items[i];
          if (mask >= 0x40) { a = 1; mask >>= 6; }
          while (!(mask & 1)) { mask >>= 1; b++; }
          tab.push_back(const_cast<Scalar*>(this->solutions[i]->get_values(this->num_components == 1 ? a : j, b)));
          if (tab[i] == nullptr)
            throw Hermes::Exceptions::Exception("Value of 'item%d' is incorrect in filter definition.", i + 1);
        }

        // apply the filter
        filter_fn(np, tab, this->values[j][0]);
      }
      this->values_valid = true;
    }

    template<typename Scalar>
    Func<Scalar>* SimpleFilter<Scalar>::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      std::vector<Scalar> val;
      for (int i = 0; i < this->solutions.size(); i++)
        val.push_back(this->solutions[i]->get_pt_value(x, y, use_MeshHashGrid, e)->val[0]);

      Func<Scalar>* toReturn = new Func<Scalar>(1, 1);

      Scalar result;
      std::vector<const Scalar*> values;
      for (int i = 0; i < this->solutions.size(); i++)
        values.push_back(&val[i]);

      // apply the filter
      filter_fn(1, values, &result);

      toReturn->val[0] = result;
      return toReturn;
    }

    ComplexFilter::ComplexFilter() : Filter<double>()
    {
      this->unimesh = false;
    }

    ComplexFilter::ComplexFilter(MeshFunctionSharedPtr<std::complex<double> > solution, int item) : Filter<double>()
    {
      this->unimesh = false;
      this->sln_complex = solution;
      this->num_components = solution->get_num_components();
      this->mesh = solution->get_mesh();
      this->solutions_sub_idx.push_back(0);
      set_quad_2d(&g_quad_2d_std);
    }

    ComplexFilter::~ComplexFilter()
    {
      this->free();
    }

    void ComplexFilter::free()
    {
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
      this->solutions_sub_idx[0] = 0;

      // fixme
      this->order = 20;
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

    void ComplexFilter::precalculate(unsigned short order, unsigned short mask)
    {
#ifdef H2D_USE_SECOND_DERIVATIVES
      if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
#else
      if (mask & (H2D_FN_DX | H2D_FN_DY))
#endif
        throw Hermes::Exceptions::Exception("Filter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      unsigned char np = quad->get_num_points(order, this->element->get_mode());

      this->sln_complex->set_quad_order(order, H2D_FN_VAL);

      // obtain corresponding tables
      filter_fn(np, const_cast<std::complex<double>*>(this->sln_complex->get_values(0, 0)), this->values[0][0]);
      if (num_components > 1)
        filter_fn(np, const_cast<std::complex<double>*>(this->sln_complex->get_values(1, 0)), this->values[1][0]);
      this->values_valid = true;
    }

    Func<double>* ComplexFilter::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      Func<std::complex<double> >* val = this->sln_complex->get_pt_value(x, y, use_MeshHashGrid, e);

      Func<double>* toReturn = new Func<double>(1, this->num_components);

      double result;

      // apply the filter
      filter_fn(1, val->val, &result);

      if (this->num_components == 1)
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
    DXDYFilter<Scalar>::DXDYFilter() : Filter<Scalar>()
    {
    }

    template<typename Scalar>
    DXDYFilter<Scalar>::DXDYFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions) : Filter<Scalar>(solutions)
    {
      init_components();
    }

    template<typename Scalar>
    DXDYFilter<Scalar>::~DXDYFilter()
    {
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::init_components()
    {
      this->num_components = this->solutions[0]->get_num_components();
      for (int i = 1; i < this->solutions.size(); i++)
        if (this->solutions[i]->get_num_components() != this->num_components)
          throw Hermes::Exceptions::Exception("Filter: Solutions do not have the same number of components!");
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::init(std::vector<MeshFunctionSharedPtr<Scalar> > solutions)
    {
      Filter<Scalar>::init();
      init_components();
    }

    template<typename Scalar>
    void DXDYFilter<Scalar>::precalculate(unsigned short order, unsigned short mask)
    {
      Quad2D* quad = this->quads[this->cur_quad];
      unsigned char np = quad->get_num_points(order, this->element->get_mode());

      // precalculate all solutions
      for (int i = 0; i < this->solutions.size(); i++)
        this->solutions[i]->set_quad_order(order, H2D_FN_DEFAULT);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain solution tables
        double *x, *y;
        const Scalar *val[H2D_MAX_COMPONENTS], *dx[H2D_MAX_COMPONENTS], *dy[H2D_MAX_COMPONENTS];
        x = this->solutions[0]->get_refmap()->get_phys_x(order);
        y = this->solutions[0]->get_refmap()->get_phys_y(order);

        for (int i = 0; i < this->solutions.size(); i++)
        {
          val[i] = this->solutions[i]->get_fn_values(j);
          dx[i] = this->solutions[i]->get_dx_values(j);
          dy[i] = this->solutions[i]->get_dy_values(j);
        }

        std::vector<const Scalar *> values_vector;
        std::vector<const Scalar *> dx_vector;
        std::vector<const Scalar *> dy_vector;

        for (int i = 0; i < this->solutions.size(); i++)
        {
          values_vector.push_back(val[i]);
          dx_vector.push_back(dx[i]);
          dy_vector.push_back(dy[i]);
        }

        // apply the filter
        filter_fn(np, x, y, values_vector, dx_vector, dy_vector, this->values[j][0], this->values[j][1], this->values[j][2]);
      }

      this->values_valid = true;
    }

    template<typename Scalar>
    Func<Scalar>* DXDYFilter<Scalar>::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      this->warn("DXDYFilter<Scalar>::get_pt_value not implemented.");
      return 0;
    }

    template<typename Scalar>
    void MagFilter<Scalar>::filter_fn(int n, const std::vector<const Scalar*>& values, Scalar* result)
    {
      for (int i = 0; i < n; i++)

      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          result[i] += sqr(values[j][i]);
        result[i] = sqrt(result[i]);
      }
    };

    template<typename Scalar>
    MagFilter<Scalar>::MagFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions, std::vector<int> items) : SimpleFilter<Scalar>(solutions, items)
    {
      Filter<Scalar>::init();
    };

    template<typename Scalar>
    MagFilter<Scalar>::MagFilter(MeshFunctionSharedPtr<Scalar> sln1, int item1)
      : SimpleFilter<Scalar>()
    {
      this->solutions = { sln1, sln1 };
      this->items = { item1 & H2D_FN_COMPONENT_0, item1 & H2D_FN_COMPONENT_1 };
      if (sln1->get_num_components() < 2)
        throw Hermes::Exceptions::Exception("The single-argument constructor is intended for vector-valued solutions.");
      Filter<Scalar>::init();
    };

    template<typename Scalar>
    MagFilter<Scalar>::~MagFilter()
    {
    };

    template<typename Scalar>
    MeshFunction<Scalar>* MagFilter<Scalar>::clone() const
    {
      std::vector<MeshFunctionSharedPtr<Scalar> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      MagFilter<Scalar>* filter = new MagFilter<Scalar>(slns, items);
      return filter;
    }

    void TopValFilter::filter_fn(int n, const std::vector<const double*>& values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          if (values[j][i] > limits[j])
            result[i] = limits[j];
          else
            result[i] = values[j][i];
      }
    };

    TopValFilter::TopValFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, std::vector<double> limits, std::vector<int> items) : SimpleFilter<double>(solutions, items), limits(limits)
    {
    };

    TopValFilter::TopValFilter(MeshFunctionSharedPtr<double> sln, double limit, int item)
      : SimpleFilter<double>()
    {
      this->limits.push_back(limit);
      this->solutions.push_back(sln);
      this->items.push_back(item);
      Filter<double>::init();
    };

    TopValFilter::~TopValFilter()
    {
    }

    MeshFunction<double>* TopValFilter::clone() const
    {
      std::vector<MeshFunctionSharedPtr<double> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      TopValFilter* filter = new TopValFilter(slns, limits, items);
      return filter;
    }

    void BottomValFilter::filter_fn(int n, const std::vector<const double*>& values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          if (values[j][i] < limits[j])
            result[i] = limits[j];
          else
            result[i] = values[j][i];
      }
    };

    BottomValFilter::BottomValFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, std::vector<double> limits, std::vector<int> items) : SimpleFilter<double>(solutions, items), limits(limits)
    {
    };

    BottomValFilter::BottomValFilter(MeshFunctionSharedPtr<double> sln, double limit, int item)
      : SimpleFilter<double>()
    {
      this->limits.push_back(limit);
      this->solutions.push_back(sln);
      this->items.push_back(item);
      Filter<double>::init();
    };

    BottomValFilter::~BottomValFilter()
    {
    }

    MeshFunction<double>* BottomValFilter::clone() const
    {
      std::vector<MeshFunctionSharedPtr<double> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      BottomValFilter* filter = new BottomValFilter(slns, limits, items);
      return filter;
    }

    void ValFilter::filter_fn(int n, const std::vector<const double*>& values, double* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          if (values[j][i] < low_limits[j])
            result[i] = low_limits[j];
          else
            if (values[j][i] > high_limits[j])
              result[i] = high_limits[j];
            else
              result[i] = values[j][i];
      }
    };

    ValFilter::ValFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, std::vector<double> low_limits, std::vector<double> high_limits, std::vector<int> items) : SimpleFilter<double>(solutions, items), low_limits(low_limits), high_limits(high_limits)
    {
    };

    ValFilter::ValFilter(MeshFunctionSharedPtr<double> sln, double low_limit, double high_limit, int item)
      : SimpleFilter<double>()
    {
      this->low_limits.push_back(low_limit);
      this->high_limits.push_back(high_limit);
      this->solutions.push_back(sln);
      this->items.push_back(item);
      Filter<double>::init();
    };

    ValFilter::~ValFilter()
    {
    }

    MeshFunction<double>* ValFilter::clone() const
    {
      std::vector<MeshFunctionSharedPtr<double> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      ValFilter* filter = new ValFilter(slns, low_limits, high_limits, items);
      return filter;
    }

    template<typename Scalar>
    void DiffFilter<Scalar>::filter_fn(int n, const std::vector<const Scalar*>& values, Scalar* result)
    {
      for (int i = 0; i < n; i++) result[i] = values[0][i] - values.at(1)[i];
    };

    template<typename Scalar>
    DiffFilter<Scalar>::DiffFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions, std::vector<int> items) : SimpleFilter<Scalar>(solutions, items) {}

    template<typename Scalar>
    DiffFilter<Scalar>::~DiffFilter()
    {
    }

    template<typename Scalar>
    MeshFunction<Scalar>* DiffFilter<Scalar>::clone() const
    {
      std::vector<MeshFunctionSharedPtr<Scalar> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      DiffFilter* filter = new DiffFilter<Scalar>(slns, items);
      return filter;
    }

    template<typename Scalar>
    void SumFilter<Scalar>::filter_fn(int n, const std::vector<const Scalar*>& values, Scalar* result)
    {
      for (int i = 0; i < n; i++)
      {
        result[i] = 0;
        for (unsigned int j = 0; j < values.size(); j++)
          result[i] += values[j][i];
      }
    };

    template<typename Scalar>
    SumFilter<Scalar>::SumFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions, std::vector<int> items) : SimpleFilter<Scalar>(solutions, items) {}

    template<typename Scalar>
    SumFilter<Scalar>::~SumFilter()
    {
    }

    template<typename Scalar>
    MeshFunction<Scalar>* SumFilter<Scalar>::clone() const
    {
      std::vector<MeshFunctionSharedPtr<Scalar> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      SumFilter<Scalar>* filter = new SumFilter<Scalar>(slns, items);
      return filter;
    }

    template<>
    void SquareFilter<double>::filter_fn(int n, const std::vector<const double *>& v1, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = sqr(v1[0][i]);
    };

    template<>
    void SquareFilter<std::complex<double> >::filter_fn(int n, const std::vector<const std::complex<double> *>& v1, std::complex<double> * result)
    {
      for (int i = 0; i < n; i++)
        result[i] = std::norm(v1[0][i]);
    };

    template<typename Scalar>
    SquareFilter<Scalar>::SquareFilter(std::vector<MeshFunctionSharedPtr<Scalar> > solutions, std::vector<int> items)
      : SimpleFilter<Scalar>(solutions, items)
    {
      if (solutions.size() > 1)
        throw Hermes::Exceptions::Exception("SquareFilter only supports one MeshFunction.");
    };

    template<typename Scalar>
    SquareFilter<Scalar>::~SquareFilter()
    {
    }

    template<typename Scalar>
    MeshFunction<Scalar>* SquareFilter<Scalar>::clone() const
    {
      std::vector<MeshFunctionSharedPtr<Scalar> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      SquareFilter<Scalar>* filter = new SquareFilter<Scalar>(slns, items);
      return filter;
    }

    void AbsFilter::filter_fn(int n, const std::vector<const double*>& v1, double * result)
    {
      for (int i = 0; i < n; i++)
        result[i] = std::abs(v1[0][i]);
    };

    AbsFilter::AbsFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, std::vector<int> items)
      : SimpleFilter<double>(solutions, items)
    {
      if (solutions.size() > 1)
        throw Hermes::Exceptions::Exception("AbsFilter only supports one MeshFunction.");
    };

    AbsFilter::AbsFilter(MeshFunctionSharedPtr<double> solution)
      : SimpleFilter<double>()
    {
      this->solutions.push_back(solution);

      this->items.push_back(H2D_FN_VAL);

      this->init();
      init_components();
    };

    AbsFilter::~AbsFilter()
    {
    }

    MeshFunction<double>* AbsFilter::clone() const
    {
      std::vector<MeshFunctionSharedPtr<double> > slns;
      std::vector<int> items;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        slns.push_back(this->solutions[i]->clone());
        items.push_back(this->items[i]);
      }
      AbsFilter* filter = new AbsFilter(slns, items);
      return filter;
    }

    void RealFilter::filter_fn(int n, const std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = values[i].real();
    };

    MeshFunction<double>* RealFilter::clone() const
    {
      RealFilter* filter = new RealFilter(this->sln_complex->clone(), this->item);
      return filter;
    }

    RealFilter::RealFilter()
      : ComplexFilter()
    {
    };

    RealFilter::RealFilter(MeshFunctionSharedPtr<std::complex<double> > solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    RealFilter::~RealFilter()
    {
    }

    void ImagFilter::filter_fn(int n, const std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = values[i].imag();
    };

    ImagFilter::ImagFilter(MeshFunctionSharedPtr<std::complex<double> > solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    ImagFilter::~ImagFilter()
    {
    }

    MeshFunction<double>* ImagFilter::clone() const
    {
      ImagFilter* filter = new ImagFilter(this->sln_complex->clone(), this->item);
      return filter;
    }

    void ComplexAbsFilter::filter_fn(int n, const std::complex<double>* values, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = sqrt(sqr(values[i].real()) + sqr(values[i].imag()));
    };

    MeshFunction<double>* ComplexAbsFilter::clone() const
    {
      ComplexAbsFilter* filter = new ComplexAbsFilter(this->sln_complex->clone(), this->item);
      return filter;
    }

    ComplexAbsFilter::ComplexAbsFilter(MeshFunctionSharedPtr<std::complex<double> > solution, int item)
      : ComplexFilter(solution, item)
    {
    };

    ComplexAbsFilter::~ComplexAbsFilter()
    {
    }

    void AngleFilter::filter_fn(int n, const std::vector<const std::complex<double>*>& v1, double* result)
    {
      for (int i = 0; i < n; i++)
        result[i] = atan2(v1[0][i].imag(), v1[0][i].real());
    };

    AngleFilter::AngleFilter(std::vector<MeshFunctionSharedPtr<std::complex<double> > > solutions, std::vector<int> items)
      : SimpleFilter<std::complex<double> >(solutions, items)
    {
      if (solutions.size() > 1)
        throw Hermes::Exceptions::Exception("RealFilter only supports one MeshFunction.");
    };

    AngleFilter::~AngleFilter()
    {
    }

    void VonMisesFilter::precalculate(unsigned short order, unsigned short mask)
    {
#ifdef H2D_USE_SECOND_DERIVATIVES
      if (mask & (H2D_FN_DX | H2D_FN_DY | H2D_FN_DXX | H2D_FN_DYY | H2D_FN_DXY))
#else
      if (mask & (H2D_FN_DX | H2D_FN_DY))
#endif
        throw Hermes::Exceptions::Exception("VonMisesFilter not defined for derivatives.");

      Quad2D* quad = this->quads[this->cur_quad];
      unsigned char np = quad->get_num_points(order, this->element->get_mode());

      this->solutions[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
      this->solutions[1]->set_quad_order(order, H2D_FN_DX | H2D_FN_DY);

      const double *dudx = this->solutions[0]->get_dx_values();
      const double *dudy = this->solutions[0]->get_dy_values();
      const double *dvdx = this->solutions[1]->get_dx_values();
      const double *dvdy = this->solutions[1]->get_dy_values();
      const double *uval = this->solutions[0]->get_fn_values();
      update_refmap();
      double *x = refmap.get_phys_x(order);

      for (int i = 0; i < np; i++)
      {
        // stress tensor
        double tz = lambda*(dudx[i] + dvdy[i]);
        double tx = tz + 2 * mu*dudx[i];
        double ty = tz + 2 * mu*dvdy[i];
        if (cyl) tz += 2 * mu*uval[i] / x[i];
        double txy = mu*(dudy[i] + dvdx[i]);

        // Von Mises stress
        this->values[0][0][i] = 1.0 / sqrt(2.0) * sqrt(sqr(tx - ty) + sqr(ty - tz) + sqr(tz - tx) + 6 * sqr(txy));
      }
      this->values_valid = true;
    }

    Func<double>* VonMisesFilter::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      this->warn("VonMisesFilter<Scalar>::get_pt_value not implemented.");
      return 0;
    }

    VonMisesFilter::VonMisesFilter(std::vector<MeshFunctionSharedPtr<double> > solutions, double lambda, double mu,
      int cyl, int item1, int item2)
      : Filter<double>(solutions)
    {
      this->mu = mu;
      this->lambda = lambda;
      this->cyl = cyl;
      this->item1 = item1;
      this->item2 = item2;
    }

    VonMisesFilter::~VonMisesFilter()
    {
    }

    MeshFunction<double>* VonMisesFilter::clone() const
    {
      std::vector<MeshFunctionSharedPtr<double> > slns;
      for (int i = 0; i < this->solutions.size(); i++)
        slns.push_back(this->solutions[i]->clone());
      VonMisesFilter* filter = new VonMisesFilter(slns, lambda, mu, cyl, item1, item2);
      return filter;
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::precalculate(unsigned short order, unsigned short mask)
    {
      Quad2D* quad = this->quads[this->cur_quad];
      unsigned char np = quad->get_num_points(order, this->element->get_mode());
      struct Filter<Scalar>::Node* node = this->new_node(H2D_FN_DEFAULT, np);

      // precalculate all solutions
      for (int i = 0; i < this->solutions.size(); i++)
        this->solutions[i]->set_quad_order(order);

      for (int j = 0; j < this->num_components; j++)
      {
        // obtain solution tables
        std::vector<Scalar*> val, dx, dy;
        for (int i = 0; i < this->solutions.size(); i++)
        {
          val.push_back(this->solutions[i]->get_fn_values(j));
          dx.push_back(this->solutions[i]->get_dx_values(j));
          dy.push_back(this->solutions[i]->get_dy_values(j));
        }
        if (this->solutions.size() == 2)
        {
          for (int i = 0; i < np; i++)
          {
            node->values[j][0][i] = tau_frac * (val[1][i] - val[0][i]) + val[1][i];
            node->values[j][1][i] = tau_frac * (dx[1][i] - dx[0][i]) + dx[1][i];
            node->values[j][2][i] = tau_frac * (dy[1][i] - dy[0][i]) + dy[1][i];
          }
        }
        else
        {
          for (int i = 0; i < np; i++)
          {
            node->values[j][0][i] = val[0][i];
            node->values[j][1][i] = dx[0][i];
            node->values[j][2][i] = dy[0][i];
          }
        }
      }
      this->values_valid = true;
    }

    template<typename Scalar>
    Func<Scalar>* LinearFilter<Scalar>::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      this->warn("LinearFilter<Scalar>::get_pt_value not implemented.");
      return 0;
    }

    template<typename Scalar>
    LinearFilter<Scalar>::LinearFilter(MeshFunctionSharedPtr<Scalar> old) // : Filter(old)
    {
      init_components();
    }

    template<typename Scalar>
    LinearFilter<Scalar>::LinearFilter(MeshFunctionSharedPtr<Scalar> older, MeshFunctionSharedPtr<Scalar> old, double tau_frac)
      : Filter<Scalar>({ older, old })
    {
      this->tau_frac = tau_frac;
      init_components();
    }

    template<typename Scalar>
    LinearFilter<Scalar>::~LinearFilter()
    {
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::init_components()
    {
      this->num_components = this->solutions[0]->get_num_components();
      for (int i = 1; i < this->solutions.size(); i++)
        if (this->solutions[i]->get_num_components() != this->num_components)
          throw Hermes::Exceptions::Exception("Filter: Solutions do not have the same number of components!");
    }

    template<typename Scalar>
    void LinearFilter<Scalar>::set_active_element(Element* e)
    {
      Filter<Scalar>::set_active_element(e);

      this->order = 0;
      for (int i = 0; i < this->solutions.size(); i++)
      {
        int o = this->solutions[i]->get_fn_order();
        if (o > this->order) this->order = o;
      }
    }

    template class HERMES_API Filter < double > ;
    template class HERMES_API Filter < std::complex<double> > ;
    template class HERMES_API SimpleFilter < double > ;
    template class HERMES_API SimpleFilter < std::complex<double> > ;
    template class HERMES_API DXDYFilter < double > ;
    template class HERMES_API DXDYFilter < std::complex<double> > ;
    template class HERMES_API MagFilter < double > ;
    template class HERMES_API MagFilter < std::complex<double> > ;
    template class HERMES_API DiffFilter < double > ;
    template class HERMES_API DiffFilter < std::complex<double> > ;
    template class HERMES_API SumFilter < double > ;
    template class HERMES_API SumFilter < std::complex<double> > ;
    template class HERMES_API SquareFilter < double > ;
    template class HERMES_API SquareFilter < std::complex<double> > ;
  }
}