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

#include "exact_solution.h"
#include "forms.h"
#include "solution_h2d_xml.h"
#include "ogprojection.h"
#include "api2d.h"
#include <algorithm>

namespace Hermes
{
  namespace Hermes2D
  {
    static double3* cheb_tab_tri[11];
    static double3* cheb_tab_quad[11];
    static int      cheb_np_tri[11];
    static int      cheb_np_quad[11];

    static double3** cheb_tab[2] = { cheb_tab_tri, cheb_tab_quad };
    static int*      cheb_np[2]  = { cheb_np_tri,  cheb_np_quad  };

    static class Quad2DCheb : public Quad2D
    {
    public:

      Quad2DCheb()
      {
        max_order[0]  = max_order[1]  = 10;
        num_tables[0] = num_tables[1] = 11;
        tables = cheb_tab;
        np = cheb_np;

        tables[0][0] = tables[1][0] = NULL;
        np[0][0] = np[1][0] = 0;

        int i, j, k, n, m;
        double3* pt;
        for (int mode_i = 0; mode_i <= 1; mode_i++)
        {
          for (k = 0; k <= 10; k++)
          {
            np[mode_i][k] = n = mode_i ? sqr(k + 1) : (k + 1)*(k + 2)/2;
            tables[mode_i][k] = pt = new double3[n];

            for (i = k, m = 0; i >= 0; i--)
              for (j = k; j >= (mode_i ? 0 : k-i); j--, m++) {
                pt[m][0] = k ? cos(j * M_PI / k) : 1.0;
                pt[m][1] = k ? cos(i * M_PI / k) : 1.0;
                pt[m][2] = 1.0;
              }
          }
        }
      };

      ~Quad2DCheb()
      {
        for (int mode_i = 0; mode_i <= 1; mode_i++)
          for (int k = 1; k <= 10; k++)
            delete [] tables[mode_i][k];
      }

      virtual void dummy_fn() {}
    } g_quad_2d_cheb;

    template<typename Scalar>
    void Solution<Scalar>::set_static_verbose_output(bool verbose)
    {
      Solution<Scalar>::static_verbose_output = verbose;
    }

    template<typename Scalar>
    bool Solution<Scalar>::static_verbose_output = false;

    template<typename Scalar>
    void Solution<Scalar>::init()
    {
      memset(elems,  0, sizeof(elems));
      memset(oldest, 0, sizeof(oldest));
      transform = true;
      sln_type = HERMES_UNDEF;
      this->num_components = 0;
      e_last = NULL;

      mono_coeffs = NULL;
      elem_coeffs[0] = elem_coeffs[1] = NULL;
      elem_orders = NULL;
      dxdy_buffer = NULL;
      num_coeffs = num_elems = 0;
      num_dofs = -1;

      this->set_quad_2d(&g_quad_2d_std);
    }

    template<>
    Solution<double>::Solution()
      : MeshFunction<double>()
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
    }

    template<>
    Solution<std::complex<double> >::Solution()
      : MeshFunction<std::complex<double> >()
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
    }

    template<>
    Solution<double>::Solution(MeshSharedPtr mesh) : MeshFunction<double>(mesh)
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
      this->mesh = mesh;
    }

    template<>
    Solution<std::complex<double> >::Solution(MeshSharedPtr mesh) : MeshFunction<std::complex<double> >(mesh)
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
      this->mesh = mesh;
    }

    template<>
    Solution<double>::Solution(SpaceSharedPtr<double> s, Vector<double>* coeff_vec) : MeshFunction<double>(s->get_mesh())
    {
      space_type = s->get_type();
      this->init();
      this->mesh = s->get_mesh();
      Solution<double>::vector_to_solution(coeff_vec, s, this);
    }

    template<>
    Solution<std::complex<double> >::Solution(SpaceSharedPtr<std::complex<double> > s, Vector<std::complex<double> >* coeff_vec) : MeshFunction<std::complex<double> >(s->get_mesh())
    {
      space_type = s->get_type();
      this->init();
      this->mesh = s->get_mesh();
      Solution<std::complex<double> >::vector_to_solution(coeff_vec, s, this);
    }

    template<>
    Solution<double>::Solution(SpaceSharedPtr<double> s, double* coeff_vec) : MeshFunction<double>(s->get_mesh())
    {
      space_type = s->get_type();
      this->init();
      this->mesh = s->get_mesh();
      Solution<double>::vector_to_solution(coeff_vec, s, this);
    }

    template<>
    Solution<std::complex<double> >::Solution(SpaceSharedPtr<std::complex<double> > s, std::complex<double> * coeff_vec) : MeshFunction<std::complex<double> >(s->get_mesh())
    {
      space_type = s->get_type();
      this->init();
      this->mesh = s->get_mesh();
      Solution<std::complex<double> >::vector_to_solution(coeff_vec, s, this);
    }

    template<typename Scalar>
    void Solution<Scalar>::copy(const MeshFunction<Scalar>* sln)
    {
      const Solution<Scalar>* solution = dynamic_cast<const Solution<Scalar>*>(sln);
      if(solution == NULL)
        throw Exceptions::Exception("The instance is in fact not a Solution instance in copy().");

      if(solution->sln_type == HERMES_UNDEF) 
        throw Hermes::Exceptions::Exception("Solution being copied is uninitialized.");
      free();

      this->mesh = solution->mesh;

      sln_type = solution->sln_type;
      space_type = solution->get_space_type();
      this->num_components = solution->num_components;
      num_dofs = solution->num_dofs;

      if(solution->sln_type == HERMES_SLN) // standard solution: copy coefficient arrays
      {
        num_coeffs = solution->num_coeffs;
        num_elems = solution->num_elems;

        mono_coeffs = new Scalar[num_coeffs];
        memcpy(mono_coeffs, solution->mono_coeffs, sizeof(Scalar) * num_coeffs);

        for (int l = 0; l < this->num_components; l++)
        {
          elem_coeffs[l] = new int[num_elems];
          memcpy(elem_coeffs[l], solution->elem_coeffs[l], sizeof(int) * num_elems);
        }

        elem_orders = new int[num_elems];
        memcpy(elem_orders, solution->elem_orders, sizeof(int) * num_elems);

        init_dxdy_buffer();
      }
      else // Const, exact handled differently.
        throw Hermes::Exceptions::Exception("Undefined or exact solutions cannot be copied into an instance of Solution already coming from computation.");

      this->element = NULL;
    }

    template<typename Scalar>
    MeshFunction<Scalar>* Solution<Scalar>::clone() const
    {
      Solution<Scalar>* sln = new Solution<Scalar>();
      sln->copy(this);
      return sln;
    }

    template<typename Scalar>
    void Solution<Scalar>::free_tables()
    {
      for (int i = 0; i < H2D_MAX_QUADRATURES; i++)
        for (int j = 0; j < H2D_SOLUTION_ELEMENT_CACHE_SIZE; j++)
        {
          tables[i][j].run_for_all(Function<Scalar>::Node::DeallocationFunction);
          tables[i][j].clear();
          elems[i][j] = NULL;
        }
    }

    template<typename Scalar>
    void Solution<Scalar>::free()
    {
      if(mono_coeffs  != NULL) { delete [] mono_coeffs;   mono_coeffs = NULL;  }
      if(elem_orders != NULL) { delete [] elem_orders;  elem_orders = NULL; }
      if(dxdy_buffer != NULL) { delete [] dxdy_buffer;  dxdy_buffer = NULL; }

      for (int i = 0; i < this->num_components; i++)
        if(elem_coeffs[i] != NULL)
        { delete [] elem_coeffs[i];  elem_coeffs[i] = NULL; }

        e_last = NULL;

        free_tables();
    }

    template<typename Scalar>
    Solution<Scalar>::~Solution()
    {
      free();
      space_type = HERMES_INVALID_SPACE;
    }

    static class mono_lu_init
    {
    public:

      // this is a set of LU-decomposed matrices shared by all Solutions
      double** mat[2][11];
      int* perm[2][11];

      mono_lu_init()
      {
        memset(mat, 0, sizeof(mat));
      }

      ~mono_lu_init()
      {
        for (int m = 0; m <= 1; m++)
          for (int i = 0; i <= 10; i++)
            if(mat[m][i] != NULL) {
              delete [] mat[m][i];
              delete [] perm[m][i];
            }
      }
    }
    mono_lu;

    template<typename Scalar>
    double** Solution<Scalar>::calc_mono_matrix(int o, int*& perm)
    {
      int i, j, k, l, m, row;
      double x, y, xn, yn;
      int n = this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;

      // loop through all chebyshev points
      double** mat = new_matrix<double>(n, n);
      for (k = o, row = 0; k >= 0; k--)
      {
        y = o ? cos(k * M_PI / o) : 1.0;
        for (l = o; l >= (this->mode ? 0 : o-k); l--, row++)
        {
          x = o ? cos(l * M_PI / o) : 1.0;

          // each row of the matrix contains all the monomials x^i*y^j
          for (i = 0, yn = 1.0, m = n-1;  i <= o;  i++, yn *= y)
            for (j = (this->mode ? 0 : i), xn = 1.0;  j <= o;  j++, xn *= x, m--)
              mat[row][m] = xn * yn;
        }
      }

      double d;
      perm = new int[n];
      ludcmp(mat, n, perm, &d);
      return mat;
    }

    template<typename Scalar>
    void Solution<Scalar>::set_coeff_vector(SpaceSharedPtr<Scalar> space, const Vector<Scalar>* vec,
      bool add_dir_lift, int start_index)
    {
      // Sanity check.
      if(vec == NULL) throw Exceptions::NullException(2);

      space_type = space->get_type();
      Scalar* coeffs = new Scalar[vec->length()];
      vec->extract(coeffs);
      this->set_coeff_vector(space, coeffs, add_dir_lift, start_index);
      delete [] coeffs;
    }

    template<typename Scalar>
    void Solution<Scalar>::set_coeff_vector(SpaceSharedPtr<Scalar> space, const Scalar* coeffs,
      bool add_dir_lift, int start_index)
    {
      // Initialize precalc shapeset using the space's shapeset.
      Shapeset *shapeset = space->shapeset;
      if(space->shapeset == NULL)
        throw Exceptions::Exception("Space->shapeset == NULL in Solution<Scalar>::set_coeff_vector().");
      PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
      if(pss == NULL) throw Exceptions::Exception("PrecalcShapeset could not be allocated in Solution<Scalar>::set_coeff_vector().");
      set_coeff_vector(space, pss, coeffs, add_dir_lift, start_index);
      delete pss;
    }

    template<typename Scalar>
    void Solution<Scalar>::set_coeff_vector(SpaceSharedPtr<Scalar> space, PrecalcShapeset* pss,
      const Scalar* coeff_vec, bool add_dir_lift, int start_index)
    {
      int o;
      if(Solution<Scalar>::static_verbose_output)
        Hermes::Mixins::Loggable::Static::info("Solution: set_coeff_vector called.");
      // Sanity checks.
      if(space->get_mesh() == NULL)
        throw Exceptions::Exception("Mesh == NULL in Solution<Scalar>::set_coeff_vector().");
      if(pss == NULL) throw Exceptions::NullException(2);
      if(coeff_vec == NULL) throw Exceptions::NullException(3);
      if(coeff_vec == NULL) throw Exceptions::Exception("Coefficient vector == NULL in Solution<Scalar>::set_coeff_vector().");
      if(!space->is_up_to_date())
        throw Exceptions::Exception("Provided 'space' is not up to date.");
      if(space->shapeset != pss->shapeset)
        throw Exceptions::Exception("Provided 'space' and 'pss' must have the same shapesets.");

      if(Solution<Scalar>::static_verbose_output)
        Hermes::Mixins::Loggable::Static::info("Solution: set_coeff_vector - solution being freed.");

      free();

      this->space_type = space->get_type();

      this->num_components = pss->get_num_components();
      sln_type = HERMES_SLN;

      // Copy the mesh.
      this->mesh = space->get_mesh();

      // Allocate the coefficient arrays.
      num_elems = this->mesh->get_max_element_id();
      if(elem_orders != NULL)
        delete [] elem_orders;
      elem_orders = new int[num_elems];
      memset(elem_orders, 0, sizeof(int) * num_elems);
      for (int l = 0; l < this->num_components; l++)
      {
        if(elem_coeffs[l] != NULL)
          delete [] elem_coeffs[l];
        elem_coeffs[l] = new int[num_elems];
        memset(elem_coeffs[l], 0, sizeof(int) * num_elems);
      }

      // Obtain element orders, allocate mono_coeffs.
      Element* e;
      num_coeffs = 0;
      for_all_active_elements(e, this->mesh)
      {
        this->mode = e->get_mode();
        o = space->get_element_order(e->id);
        o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
        for (unsigned int k = 0; k < e->get_nvert(); k++)
        {
          int eo = space->get_edge_order(e, k);
          if(eo > o) o = eo;
        }

        // Hcurl and Hdiv: actual order of functions is one higher than element order
        if((space->shapeset)->get_num_components() == 2) o++;

        num_coeffs += this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;
        elem_orders[e->id] = o;
      }
      num_coeffs *= this->num_components;
      if(mono_coeffs != NULL)
        delete [] mono_coeffs;
      mono_coeffs = new Scalar[num_coeffs];

      // Express the solution on elements as a linear combination of monomials.
      Quad2D* quad = &g_quad_2d_cheb;
      pss->set_quad_2d(quad);
      Scalar* mono = mono_coeffs;
      for_all_active_elements(e, this->mesh)
      {
        this->mode = e->get_mode();
        o = elem_orders[e->id];
        int np = quad->get_num_points(o, e->get_mode());

        AsmList<Scalar> al;
        space->get_element_assembly_list(e, &al);
        pss->set_active_element(e);

        for (int l = 0; l < this->num_components; l++)
        {
          // Obtain solution values for the current element.
          Scalar* val = mono;
          elem_coeffs[l][e->id] = (int) (mono - mono_coeffs);
          memset(val, 0, sizeof(Scalar)*np);
          for (unsigned int k = 0; k < al.cnt; k++)
          {
            pss->set_active_shape(al.idx[k]);
            pss->set_quad_order(o, H2D_FN_VAL);
            int dof = al.dof[k];
            double dir_lift_coeff = add_dir_lift ? 1.0 : 0.0;
            // By subtracting space->first_dof we make sure that it does not matter where the
            // enumeration of dofs in the space starts. This ca be either zero or there can be some
            // offset. By adding start_index we move to the desired section of coeff_vec.
            Scalar coef = al.coef[k] * (dof >= 0 ? coeff_vec[dof  - space->first_dof + start_index] : dir_lift_coeff);
            double* shape = pss->get_fn_values(l);
            for (int i = 0; i < np; i++)
              val[i] += shape[i] * coef;
          }
          mono += np;

          // solve for the monomial coefficients
          if(mono_lu.mat[this->mode][o] == NULL)
            mono_lu.mat[this->mode][o] = calc_mono_matrix(o, mono_lu.perm[this->mode][o]);
          lubksb(mono_lu.mat[this->mode][o], np, mono_lu.perm[this->mode][o], val);
        }
      }

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL");
      init_dxdy_buffer();
      this->element = NULL;
      if(Solution<Scalar>::static_verbose_output)
        Hermes::Mixins::Loggable::Static::info("Solution: set_coeff_vector - done.");
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Scalar* solution_vector,
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions,
      Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) 
        throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

      // If start indices are not given, calculate them using the dimension of each space.
      Hermes::vector<int> start_indices_new;
      if(start_indices.empty())
      {
        int counter = 0;
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(counter);
          counter += spaces[i]->get_num_dofs();
        }
      }
      else
      {
        if(start_indices.size() != spaces.size()) throw Hermes::Exceptions::Exception("Mismatched start indices in vector_to_solutions().");
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(start_indices[i]);
        }
      }

      for (int i=0; i < spaces.size(); i++)
      {
        if(Solution<Scalar>::static_verbose_output)
          Hermes::Mixins::Loggable::Static::info("Vector to Solution: %d-th solution", i);

        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Scalar* solution_vector, SpaceSharedPtr<Scalar> space,
      Solution<Scalar>* solution, bool add_dir_lift, int start_index)
    {
      // Sanity checks.
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);

      solution->set_coeff_vector(space, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Scalar* solution_vector, SpaceSharedPtr<Scalar> space,
      MeshFunctionSharedPtr<Scalar> solution, bool add_dir_lift, int start_index)
    {
      // Sanity checks.
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);

      Solution<Scalar>* sln = dynamic_cast<Solution<Scalar>*>(solution.get());
      if(sln == NULL)
        throw Exceptions::Exception("Passed solution is in fact not a Solution instance in vector_to_solution().");

      sln->set_coeff_vector(space, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Vector<Scalar>* solution_vector, Hermes::vector<SpaceSharedPtr<Scalar> > spaces,
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions, Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) 
        throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

      // If start indices are not given, calculate them using the dimension of each space.
      Hermes::vector<int> start_indices_new;
      if(start_indices.empty())
      {
        int counter = 0;
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(counter);
          counter += spaces[i]->get_num_dofs();
        }
      }
      else
      {
        if(start_indices.size() != spaces.size()) throw Hermes::Exceptions::Exception("Mismatched start indices in vector_to_solutions().");
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(start_indices[i]);
        }
      }

      for (int i=0; i < spaces.size(); i++)
      {
        if(Solution<Scalar>::static_verbose_output)
          Hermes::Mixins::Loggable::Static::info("Vector to Solution: %d-th solution", i);

        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions_common_dir_lift(const Vector<Scalar>* solution_vector, Hermes::vector<SpaceSharedPtr<Scalar> > spaces,
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions, bool add_dir_lift)
    {
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) 
        throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

      // If start indices are not given, calculate them using the dimension of each space.
      Hermes::vector<int> start_indices_new;
      int counter = 0;
      for (int i=0; i < spaces.size(); i++)
      {
        start_indices_new.push_back(counter);
        counter += spaces[i]->get_num_dofs();
      }

      for (int i=0; i < spaces.size(); i++)
      {
        if(Solution<Scalar>::static_verbose_output)
          Hermes::Mixins::Loggable::Static::info("Vector to Solution: %d-th solution", i);

        Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift, start_indices_new[i]);
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions_common_dir_lift(const Scalar* solution_vector, Hermes::vector<SpaceSharedPtr<Scalar> > spaces,
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions, bool add_dir_lift)
    {
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) 
        throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

      // If start indices are not given, calculate them using the dimension of each space.
      Hermes::vector<int> start_indices_new;
      int counter = 0;
      for (int i=0; i < spaces.size(); i++)
      {
        start_indices_new.push_back(counter);
        counter += spaces[i]->get_num_dofs();
      }

      for (int i=0; i < spaces.size(); i++)
      {
        if(Solution<Scalar>::static_verbose_output)
          Hermes::Mixins::Loggable::Static::info("Vector to Solution: %d-th solution", i);

        Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift, start_indices_new[i]);
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Vector<Scalar>* solution_vector, SpaceSharedPtr<Scalar> space,
      MeshFunctionSharedPtr<Scalar> solution, bool add_dir_lift, int start_index)
    {
      // Sanity checks.
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);

      Solution<Scalar>* sln = dynamic_cast<Solution<Scalar>*>(solution.get());
      if(sln == NULL)
        throw Exceptions::Exception("Passed solution is in fact not a Solution instance in vector_to_solution().");

      sln->set_coeff_vector(space, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Scalar* solution_vector, Hermes::vector<SpaceSharedPtr<Scalar> > spaces,
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions, Hermes::vector<PrecalcShapeset *> pss,
      Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector==NULL) 
        throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) 
        throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

      // If start indices are not given, calculate them using the dimension of each space.
      Hermes::vector<int> start_indices_new;
      if(start_indices.empty())
      {
        int counter = 0;
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(counter);
          counter += spaces[i]->get_num_dofs();
        }
      }
      else
      {
        if(start_indices.size() != spaces.size()) 
          throw Hermes::Exceptions::Exception("Mismatched start indices in vector_to_solutions().");
        for (int i=0; i < spaces.size(); i++)
        {
          start_indices_new.push_back(start_indices[i]);
        }
      }
      
      for (int i=0; i < spaces.size(); i++)
      {
        if(Solution<Scalar>::static_verbose_output)
          Hermes::Mixins::Loggable::Static::info("Vector to Solution: %d-th solution", i);

        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Scalar* solution_vector, SpaceSharedPtr<Scalar> space, MeshFunctionSharedPtr<Scalar> solution,
      PrecalcShapeset* pss, bool add_dir_lift, int start_index)
    {
      if(solution_vector == NULL) 
        throw Exceptions::NullException(1);
      if(pss == NULL) 
        throw Exceptions::NullException(4);

      Solution<Scalar>* sln = dynamic_cast<Solution<Scalar>*>(solution.get());
      if(sln == NULL)
        throw Exceptions::Exception("Passed solution is in fact not a Solution instance in vector_to_solution().");

      sln->set_coeff_vector(space, pss, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::set_dirichlet_lift(SpaceSharedPtr<Scalar> space, PrecalcShapeset* pss)
    {
      space_type = space->get_type();
      int ndof = space->get_num_dofs();
      Scalar *temp = new Scalar[ndof];
      memset(temp, 0, sizeof(Scalar)*ndof);
      bool add_dir_lift = true;
      int start_index = 0;
      this->set_coeff_vector(space, pss, temp, add_dir_lift, start_index);
      delete [] temp;
    }

    template<typename Scalar>
    void Solution<Scalar>::enable_transform(bool enable)
    {
      if(transform != enable)
        free_tables();
      transform = enable;
    }

    template<typename Scalar>
    void Solution<Scalar>::multiply(Scalar coef)
    {
      if(sln_type == HERMES_SLN)
      {
        for (int i = 0; i < num_coeffs; i++)
          mono_coeffs[i] *= coef;
      }
      else if(sln_type == HERMES_EXACT)
        dynamic_cast<ExactSolution<Scalar>* >(this)->exact_multiplicator *= coef;
      else
        throw Hermes::Exceptions::Exception("Uninitialized solution.");
    }

    template<typename Scalar>
    static void make_dx_coeffs(int mode, int o, Scalar* mono, Scalar* result)
    {
      int i, j, k;
      for (i = 0; i <= o; i++) {
        *result++= 0.0;
        k = mode ? o : i;
        for (j = 0; j < k; j++)
          *result++= (Scalar) (k-j) * mono[j];
        mono += k + 1;
      }
    }

    template<typename Scalar>
    static void make_dy_coeffs(int mode, int o, Scalar* mono, Scalar* result)
    {
      int i, j;
      if(mode) {
        for (j = 0; j <= o; j++)
          *result++= 0.0;
        for (i = 0; i < o; i++)
          for (j = 0; j <= o; j++)
            *result++= (Scalar) (o-i) * (*mono++);
      }
      else {
        for (i = 0; i <= o; i++) {
          *result++= 0.0;
          for (j = 0; j < i; j++)
            *result++= (Scalar) (o + 1-i) * (*mono++);
        }
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::init_dxdy_buffer()
    {
      if(dxdy_buffer != NULL)
      {
        delete [] dxdy_buffer;
        dxdy_buffer = NULL;
      }
      dxdy_buffer = new Scalar[this->num_components * 5 * 121];
    }

    template<typename Scalar>
    void Solution<Scalar>::set_active_element(Element* e)
    {
      if(e == this->element)
        return;

      if(!e->active) 
        throw Hermes::Exceptions::Exception("Cannot select inactive element. Wrong mesh?");

      MeshFunction<Scalar>::set_active_element(e);

      // try finding an existing table for e
      for (cur_elem = 0; cur_elem < H2D_SOLUTION_ELEMENT_CACHE_SIZE; cur_elem++)
        if(elems[this->cur_quad][cur_elem] == e)
          break;

      // if not found, free the oldest one and use its slot
      if(cur_elem >= H2D_SOLUTION_ELEMENT_CACHE_SIZE)
      {
        tables[this->cur_quad][oldest[this->cur_quad]].run_for_all(Function<Scalar>::Node::DeallocationFunction);
        tables[this->cur_quad][oldest[this->cur_quad]].clear();
        elems[this->cur_quad][oldest[this->cur_quad]] = NULL;

        cur_elem = oldest[this->cur_quad];
        if(++oldest[this->cur_quad] >= H2D_SOLUTION_ELEMENT_CACHE_SIZE)
          oldest[this->cur_quad] = 0;

        elems[this->cur_quad][cur_elem] = e;
      }

      if(sln_type == HERMES_SLN)
      {
        int o = this->order = elem_orders[this->element->id];
        int n = this->mode ? sqr(o + 1) : (o + 1)*(o + 2)/2;

        for (int i = 0, m = 0; i < this->num_components; i++)
        {
          Scalar* mono = mono_coeffs + elem_coeffs[i][e->id];
          dxdy_coeffs[i][0] = mono;

          make_dx_coeffs(this->mode, o, mono, dxdy_coeffs[i][1] = dxdy_buffer + m);  m += n;
          make_dy_coeffs(this->mode, o, mono, dxdy_coeffs[i][2] = dxdy_buffer + m);  m += n;
          make_dx_coeffs(this->mode, o, dxdy_coeffs[i][1], dxdy_coeffs[i][3] = dxdy_buffer + m);  m += n;
          make_dy_coeffs(this->mode, o, dxdy_coeffs[i][2], dxdy_coeffs[i][4] = dxdy_buffer + m);  m += n;
          make_dx_coeffs(this->mode, o, dxdy_coeffs[i][2], dxdy_coeffs[i][5] = dxdy_buffer + m);  m += n;
        }
      }
      else if(sln_type == HERMES_EXACT)
      {
        this->order = Hermes::Hermes2D::g_max_quad;
      }
      else
        throw Hermes::Exceptions::Exception("Uninitialized solution.");

      this->sub_tables = &tables[this->cur_quad][cur_elem];

      this->update_nodes_ptr();
    }

    template<typename Scalar>
    static inline void set_vec_num(int n, Scalar* y, Scalar num)
    {
      for (int i = 0; i < n; i++)
        y[i] = num;
    }

    template<typename Scalar>
    static inline void vec_x_vec_p_num(int n, Scalar* y, Scalar* x, Scalar num)
    {
      for (int i = 0; i < n; i++)
        y[i] = y[i]*x[i] + num;
    }

    template<typename Scalar>
    static inline void vec_x_vec_p_vec(int n, Scalar* y, Scalar* x, Scalar* z)
    {
      for (int i = 0; i < n; i++)
        y[i] = y[i]*x[i] + z[i];
    }

    static const int H2D_GRAD = H2D_FN_DX_0 | H2D_FN_DY_0;
    static const int H2D_SECOND = H2D_FN_DXX_0 | H2D_FN_DXY_0 | H2D_FN_DYY_0;
    static const int H2D_CURL = H2D_FN_DX | H2D_FN_DY;

    template<typename Scalar>
    void Solution<Scalar>::transform_values(int order, struct Function<Scalar>::Node* node, int newmask, int oldmask, int np)
    {
      double2x2 *mat, *m;
      int i, mstep = 0;

      // H1 space, L2 space
      if((space_type == HERMES_H1_SPACE)||(space_type == HERMES_L2_SPACE))
      {
#ifdef H2D_USE_SECOND_DERIVATIVES
        if(((newmask & H2D_SECOND) == H2D_SECOND && (oldmask & H2D_SECOND) != H2D_SECOND))
        {
          this->update_refmap();
          mat = this->refmap->get_inv_ref_map(order);
          double3x2 *mm, *mat2 = this->refmap->get_second_ref_map(order);
          for (i = 0, m = mat, mm = mat2; i < np; i++, m++, mm++)
          {
            Scalar vx = node->values[0][1][i];
            Scalar vy = node->values[0][2][i];
            Scalar vxx = node->values[0][3][i];
            Scalar vyy = node->values[0][4][i];
            Scalar vxy = node->values[0][5][i];

            node->values[0][3][i] = sqr((*m)[0][0])*vxx + 2*(*m)[0][1]*(*m)[0][0]*vxy + sqr((*m)[0][1])*vyy + (*mm)[0][0]*vx + (*mm)[0][1]*vy;   // dxx
            node->values[0][4][i] = sqr((*m)[1][0])*vxx + 2*(*m)[1][1]*(*m)[1][0]*vxy + sqr((*m)[1][1])*vyy + (*mm)[2][0]*vx + (*mm)[2][1]*vy;   // dyy
            node->values[0][5][i] = (*m)[0][0]*(*m)[1][0]*vxx + ((*m)[0][0]*(*m)[1][1] + (*m)[1][0]*(*m)[0][1])*vxy + (*m)[0][1]*(*m)[1][1]*vyy + (*mm)[1][0]*vx + (*mm)[1][1]*vy;   //dxy
          }
        }
#endif
        if((newmask & H2D_GRAD) == H2D_GRAD && (oldmask & H2D_GRAD) != H2D_GRAD)
        {
          this->update_refmap();
          mat = this->refmap->get_const_inv_ref_map();
          if(!this->refmap->is_jacobian_const()) { mat = this->refmap->get_inv_ref_map(order); mstep = 1; }

          for (i = 0, m = mat; i < np; i++, m += mstep)
          {
            Scalar vx = node->values[0][1][i];
            Scalar vy = node->values[0][2][i];
            node->values[0][1][i] = (*m)[0][0]*vx + (*m)[0][1]*vy;
            node->values[0][2][i] = (*m)[1][0]*vx + (*m)[1][1]*vy;
          }
        }
      }

      // Hcurl space
      else if(space_type == HERMES_HCURL_SPACE)
      {
        bool trans_val = false, trans_curl = false;
        if((newmask & H2D_FN_VAL) == H2D_FN_VAL && (oldmask & H2D_FN_VAL) != H2D_FN_VAL) trans_val  = true;
        if((newmask &   H2D_CURL) ==   H2D_CURL && (oldmask &   H2D_CURL) !=   H2D_CURL) trans_curl = true;

        if(trans_val || trans_curl)
        {
          this->update_refmap();
          mat = this->refmap->get_const_inv_ref_map();
          if(!this->refmap->is_jacobian_const()) { mat = this->refmap->get_inv_ref_map(order); mstep = 1; }

          for (i = 0, m = mat; i < np; i++, m += mstep)
          {
            if(trans_val)
            {
              Scalar vx = node->values[0][0][i];
              Scalar vy = node->values[1][0][i];
              node->values[0][0][i] = (*m)[0][0]*vx + (*m)[0][1]*vy;
              node->values[1][0][i] = (*m)[1][0]*vx + (*m)[1][1]*vy;
            }
            if(trans_curl)
            {
              Scalar e0x = node->values[0][1][i], e0y = node->values[0][2][i];
              Scalar e1x = node->values[1][1][i], e1y = node->values[1][2][i];
              node->values[1][1][i] = (*m)[0][0]*((*m)[1][0]*e0x + (*m)[1][1]*e1x) + (*m)[0][1]*((*m)[1][0]*e0y + (*m)[1][1]*e1y);
              node->values[0][2][i] = (*m)[1][0]*((*m)[0][0]*e0x + (*m)[0][1]*e1x) + (*m)[1][1]*((*m)[0][0]*e0y + (*m)[0][1]*e1y);
            }
          }
        }
      }

      // Hdiv space
      else if(space_type == HERMES_HDIV_SPACE)
      {
        if((newmask & H2D_FN_VAL) == H2D_FN_VAL && (oldmask & H2D_FN_VAL) != H2D_FN_VAL)
        {
          this->update_refmap();
          mat = this->refmap->get_const_inv_ref_map();
          if(!this->refmap->is_jacobian_const()) { mat = this->refmap->get_inv_ref_map(order); mstep = 1; }

          for (i = 0, m = mat; i < np; i++, m += mstep)
          {
            Scalar vx = node->values[0][0][i];
            Scalar vy = node->values[1][0][i];
            node->values[0][0][i] =   (*m)[1][1]*vx - (*m)[1][0]*vy;
            node->values[1][0][i] = - (*m)[0][1]*vx + (*m)[0][0]*vy;
          }
        }
      }
    }

    template<typename Scalar>
    void Solution<Scalar>::precalculate(int order, int mask)
    {
      int i, j, k, l;
      struct Function<Scalar>::Node* node = NULL;
      Quad2D* quad = this->quads[this->cur_quad];
      int np = quad->get_num_points(order, this->mode);

      if(sln_type == HERMES_SLN)
      {
        // if we are required to transform vectors, we must precalculate both their components
        const int H2D_GRAD = H2D_FN_DX_0 | H2D_FN_DY_0;
        const int H2D_SECOND = H2D_FN_DXX_0 | H2D_FN_DXY_0 | H2D_FN_DYY_0;
        const int H2D_CURL = H2D_FN_DX | H2D_FN_DY; // sic
        if(transform)
        {
          if(this->num_components == 1)                                            // H1 space or L2 space
          {
            if((mask & H2D_FN_DX_0)  || (mask & H2D_FN_DY_0))  mask |= H2D_GRAD;
            if((mask & H2D_FN_DXX_0)  || (mask & H2D_FN_DXY_0) || (mask & H2D_FN_DYY_0))  mask |= H2D_SECOND;
          }
          else if(space_type == HERMES_HCURL_SPACE)                                           // Hcurl space
          { if((mask & H2D_FN_VAL_0) || (mask & H2D_FN_VAL_1)) mask |= H2D_FN_VAL;
          if((mask & H2D_FN_DX_1)  || (mask & H2D_FN_DY_0))  mask |= H2D_CURL; }
          else                                                                // Hdiv space
          { if((mask & H2D_FN_VAL_0) || (mask & H2D_FN_VAL_1)) mask |= H2D_FN_VAL; }
        }

        int oldmask = (this->cur_node != NULL) ? this->cur_node->mask : 0;
        int newmask = mask | oldmask;
        node = this->new_node(newmask, np);

        // transform integration points by the current matrix
        Scalar* x = new Scalar[np];
        Scalar* y = new Scalar[np];
        Scalar* tx = new Scalar[np];
        double3* pt = quad->get_points(order, this->element->get_mode());
        for (i = 0; i < np; i++)
        {
          x[i] = pt[i][0] * this->ctm->m[0] + this->ctm->t[0];
          y[i] = pt[i][1] * this->ctm->m[1] + this->ctm->t[1];
        }

        // obtain the solution values, this is the core of the whole module
        int o = elem_orders[this->element->id];
        for (l = 0; l < this->num_components; l++)
        {
          for (k = 0; k < 6; k++)
          {
            if(newmask & this->idx2mask[k][l])
            {
              Scalar* result = node->values[l][k];
              if(oldmask & this->idx2mask[k][l])
              {
                // copy the old table if we have it already
                memcpy(result, this->cur_node->values[l][k], np * sizeof(Scalar));
              }
              else
              {
                // calculate the solution values using Horner's scheme
                Scalar* mono = dxdy_coeffs[l][k];
                for (i = 0; i <= o; i++)
                {
                  set_vec_num(np, tx, *mono++);
                  for (j = 1; j <= (this->mode ? o : i); j++)
                    vec_x_vec_p_num(np, tx, x, *mono++);

                  if(!i)
                    memcpy(result, tx, sizeof(Scalar)*np);
                  else
                    vec_x_vec_p_vec(np, result, y, tx);
                }
              }
            }
          }
        }

        delete [] x;
        delete [] y;
        delete [] tx;

        // transform gradient or vector solution, if required
        if(transform)
          transform_values(order, node, newmask, oldmask, np);
      }
      else if(sln_type == HERMES_EXACT)
      {
        if(mask & ~H2D_FN_DEFAULT)
          throw Hermes::Exceptions::Exception("Cannot obtain second derivatives of an exact solution.");
        node = this->new_node(mask = H2D_FN_DEFAULT, np);

        this->update_refmap();
        double* x = this->refmap->get_phys_x(order);
        double* y = this->refmap->get_phys_y(order);

        // evaluate the exact solution
        if(this->num_components == 1)
        {
          // untransform values
          if(!transform)
          {
            double2x2 *mat, *m;
            int mstep = 0;
            mat = this->refmap->get_const_inv_ref_map();
            if(!this->refmap->is_jacobian_const()) { mat = this->refmap->get_inv_ref_map(order); mstep = 1; }

            for (i = 0, m = mat; i < np; i++, m += mstep)
            {
              double jac = (*m)[0][0] *  (*m)[1][1] - (*m)[1][0] *  (*m)[0][1];
              Scalar val, dx = 0.0, dy = 0.0;
              val = (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_function(x[i], y[i], dx, dy);
              node->values[0][0][i] = val * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[0][1][i] = (  (*m)[1][1]*dx - (*m)[0][1]*dy) / jac * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[0][2][i] = (- (*m)[1][0]*dx + (*m)[0][0]*dy) / jac * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
            }
          }
          else
          {
            for (i = 0; i < np; i++)
            {
              Scalar val, dx = 0.0, dy = 0.0;
              val = (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_function(x[i], y[i], dx, dy);
              node->values[0][0][i] = val * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[0][1][i] = dx * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[0][2][i] = dy * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
            }
          }
        }
        else
        {
          for (i = 0; i < np; i++)
          {
            Scalar2<Scalar> dx (0.0, 0.0 ), dy ( 0.0, 0.0 );
            Scalar2<Scalar> val = (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_function(x[i], y[i], dx, dy);
            for (j = 0; j < 2; j++)
            {
              node->values[j][0][i] = val[j] * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[j][1][i] = dx[j] * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
              node->values[j][2][i] = dy[j] * (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_multiplicator;
            }
          }
        }
      }
      else
      {
        throw Hermes::Exceptions::Exception("Cannot obtain values -- uninitialized solution. The solution was either "
          "not calculated yet or you used the assignment operator which destroys "
          "the solution on its right-hand side.");
      }

      if(this->nodes->present(order))
      {
        assert(this->nodes->get(order) == this->cur_node);
        ::free(this->nodes->get(order));
      }
      this->nodes->add(node, order);
      this->cur_node = node;
    }

    template<>
    void Solution<double>::save(const char* filename) const
    {
      if(sln_type == HERMES_UNDEF)
        throw Exceptions::Exception("Cannot save -- uninitialized solution.");

      try
      {
        XMLSolution::solution xmlsolution(this->num_components, this->num_elems, this->num_coeffs, 0, 0);
        switch(this->get_space_type())
        {
        case HERMES_H1_SPACE:
          xmlsolution.space().set("h1");
          break;
        case HERMES_HCURL_SPACE:
          xmlsolution.space().set("hcurl");
          break;
        case HERMES_HDIV_SPACE:
          xmlsolution.space().set("hdiv");
          break;
        case HERMES_L2_SPACE:
          xmlsolution.space().set("l2");
          break;
        default:
          throw Exceptions::Exception("This type of solution can not be saved.");
        }

        for(unsigned int coeffs_i = 0; coeffs_i < this->num_coeffs; coeffs_i++)
          xmlsolution.mono_coeffs().push_back(XMLSolution::mono_coeffs(coeffs_i, mono_coeffs[coeffs_i]));

        for(unsigned int elems_i = 0; elems_i < this->num_elems; elems_i++)
          xmlsolution.elem_orders().push_back(XMLSolution::elem_orders(elems_i, elem_orders[elems_i]));

        for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
        {
          xmlsolution.component().push_back(XMLSolution::component());
          xmlsolution.component().back().component_number() = component_i;
          for(unsigned int elems_i = 0; elems_i < this->num_elems; elems_i++)
            xmlsolution.component().back().elem_coeffs().push_back(XMLSolution::elem_coeffs(elems_i, elem_coeffs[component_i][elems_i]));
        }

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);

        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;

        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void Solution<std::complex<double> >::save(const char* filename) const
    {
      if(sln_type == HERMES_UNDEF)
        throw Exceptions::Exception("Cannot save -- uninitialized solution.");

      try
      {
        XMLSolution::solution xmlsolution(this->num_components, this->num_elems, this->num_coeffs, 0, 1);

        switch(this->get_space_type())
        {
        case HERMES_H1_SPACE:
          xmlsolution.space().set("h1");
          break;
        case HERMES_HCURL_SPACE:
          xmlsolution.space().set("hcurl");
          break;
        case HERMES_HDIV_SPACE:
          xmlsolution.space().set("hdiv");
          break;
        case HERMES_L2_SPACE:
          xmlsolution.space().set("l2");
          break;
        default:
          throw Exceptions::Exception("This type of solution can not be saved.");
        }

        for(unsigned int coeffs_i = 0; coeffs_i < this->num_coeffs; coeffs_i++)
        {
          xmlsolution.mono_coeffs().push_back(XMLSolution::mono_coeffs(coeffs_i, mono_coeffs[coeffs_i].real()));
          xmlsolution.mono_coeffs().back().im() = mono_coeffs[coeffs_i].imag();
        }

        for(unsigned int elems_i = 0; elems_i < this->num_elems; elems_i++)
          xmlsolution.elem_orders().push_back(XMLSolution::elem_orders(elems_i, elem_orders[elems_i]));

        for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
        {
          xmlsolution.component().push_back(XMLSolution::component());
          xmlsolution.component().back().component_number() = component_i;
          for(unsigned int elems_i = 0; elems_i < this->num_elems; elems_i++)
            xmlsolution.component().back().elem_coeffs().push_back(XMLSolution::elem_coeffs(elems_i, elem_coeffs[component_i][elems_i]));
        }

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);

        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;

        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void Solution<double>::load(const char* filename, SpaceSharedPtr<double> space)
    {
      free();
      this->mesh = space->get_mesh();
      this->space_type = space->get_type();

      try
      {
        ::xml_schema::flags parsing_flags = 0;
        if(!this->validate)
          parsing_flags = xml_schema::flags::dont_validate;

        std::auto_ptr<XMLSolution::solution> parsed_xml_solution(XMLSolution::solution_(filename, parsing_flags));
        sln_type = parsed_xml_solution->exact() == 0 ? HERMES_SLN : HERMES_EXACT;

        if(sln_type == HERMES_EXACT)
        {
          switch(parsed_xml_solution->ncmp())
          {
          case 1:
            if(parsed_xml_solution->exactC() == 0)
            {
              double* coeff_vec = new double[space->get_num_dofs()];
              MeshFunctionSharedPtr<double> sln(new ConstantSolution<double>(this->mesh, parsed_xml_solution->exactCXR().get()));
              OGProjection<double>::project_global(space, sln, coeff_vec);
              this->set_coeff_vector(space, coeff_vec, true, 0);
              sln_type = HERMES_SLN;
            }
            else
              throw Hermes::Exceptions::SolutionLoadFailureException("Mismatched real - complex exact solutions.");
            break;
          case 2:
            if(parsed_xml_solution->exactC() == 0)
            {
              double* coeff_vec = new double[space->get_num_dofs()];
              MeshFunctionSharedPtr<double> sln(new ConstantSolutionVector<double>(this->mesh, parsed_xml_solution->exactCXR().get(), parsed_xml_solution->exactCYR().get()));
              OGProjection<double>::project_global(space, sln, coeff_vec);
              this->set_coeff_vector(space, coeff_vec, true, 0);
              sln_type = HERMES_SLN;
            }
            else
              throw Hermes::Exceptions::SolutionLoadFailureException("Mismatched real - complex exact solutions.");
            break;
          }
        }
        else
        {
          if(!strcmp(parsed_xml_solution->space().get().c_str(),"h1"))
            if(this->space_type != HERMES_H1_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"l2"))
            if(this->space_type != HERMES_L2_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"hcurl"))
            if(this->space_type != HERMES_HCURL_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"hdiv"))
            if(this->space_type != HERMES_HDIV_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          this->num_coeffs = parsed_xml_solution->nc();
          this->num_elems = parsed_xml_solution->nel();
          this->num_components = parsed_xml_solution->ncmp();

          this->mono_coeffs = new double[num_coeffs];
          memset(this->mono_coeffs, 0, this->num_coeffs*sizeof(double));

          for(unsigned int component_i = 0; component_i < num_components; component_i++)
            elem_coeffs[component_i] = new int[num_elems];

          this->elem_orders = new int[num_elems];

          for (unsigned int coeffs_i = 0; coeffs_i < num_coeffs; coeffs_i++)
            this->mono_coeffs[parsed_xml_solution->mono_coeffs().at(coeffs_i).id()] = parsed_xml_solution->mono_coeffs().at(coeffs_i).re();

          for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
            this->elem_orders[parsed_xml_solution->elem_orders().at(elems_i).id()] = parsed_xml_solution->elem_orders().at(elems_i).ord();

          for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
            for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
              this->elem_coeffs[component_i][parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).id()] = parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).c();

        }
        init_dxdy_buffer();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionLoadFailureException(e.what());
      }
    }

    template<>
    void Solution<std::complex<double> >::load(const char* filename, SpaceSharedPtr<std::complex<double> > space)
    {
      free();
      sln_type = HERMES_SLN;
      this->mesh = space->get_mesh();
      this->space_type = space->get_type();

      try
      {
        ::xml_schema::flags parsing_flags = 0;
        if(!this->validate)
          parsing_flags = xml_schema::flags::dont_validate;

        std::auto_ptr<XMLSolution::solution> parsed_xml_solution(XMLSolution::solution_(filename, parsing_flags));
        sln_type = parsed_xml_solution->exact() == 0 ? HERMES_SLN : HERMES_EXACT;

        if(sln_type == HERMES_EXACT)
        {
          switch(parsed_xml_solution->ncmp())
          {
          case 1:
            if(parsed_xml_solution->exactC() == 1)
            {
              std::complex<double>* coeff_vec = new std::complex<double>[space->get_num_dofs()];
              MeshFunctionSharedPtr<std::complex<double> > sln(new ConstantSolution<std::complex<double> >(this->mesh, parsed_xml_solution->exactCXR().get()));
              OGProjection<std::complex<double> >::project_global(space, sln, coeff_vec);
              this->set_coeff_vector(space, coeff_vec, true, 0);
              sln_type = HERMES_SLN;
            }
            else
              throw Hermes::Exceptions::SolutionLoadFailureException("Mismatched real - complex exact solutions.");
            break;
          case 2:
            if(parsed_xml_solution->exactC() == 1)
            {
              std::complex<double>* coeff_vec = new std::complex<double>[space->get_num_dofs()];
              MeshFunctionSharedPtr<std::complex<double> > sln(new ConstantSolutionVector<std::complex<double> >(this->mesh, std::complex<double>(parsed_xml_solution->exactCXR().get(), parsed_xml_solution->exactCXC().get()), std::complex<double>(parsed_xml_solution->exactCYR().get(), parsed_xml_solution->exactCYC().get())));
              OGProjection<std::complex<double> >::project_global(space, sln, coeff_vec);
              this->set_coeff_vector(space, coeff_vec, true, 0);
              sln_type = HERMES_SLN;
            }
            else
              throw Hermes::Exceptions::SolutionLoadFailureException("Mismatched real - complex exact solutions.");
            break;
          }
        }
        else
        {
          if(!strcmp(parsed_xml_solution->space().get().c_str(),"h1"))
            if(this->space_type != HERMES_H1_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"l2"))
            if(this->space_type != HERMES_L2_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"hcurl"))
            if(this->space_type != HERMES_HCURL_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          if(!strcmp(parsed_xml_solution->space().get().c_str(),"hdiv"))
            if(this->space_type != HERMES_HDIV_SPACE)
              throw Exceptions::Exception("Space types not compliant in Solution::load().");

          ::xml_schema::flags parsing_flags = 0;
          if(!this->validate)
            parsing_flags = xml_schema::flags::dont_validate;

          std::auto_ptr<XMLSolution::solution> parsed_xml_solution(XMLSolution::solution_(filename, parsing_flags));

          this->num_coeffs = parsed_xml_solution->nc();
          this->num_elems = parsed_xml_solution->nel();
          this->num_components = parsed_xml_solution->ncmp();

          this->mono_coeffs = new std::complex<double>[num_coeffs];
          memset(this->mono_coeffs, 0, this->num_coeffs*sizeof(std::complex<double>));

          for(unsigned int component_i = 0; component_i < num_components; component_i++)
            elem_coeffs[component_i] = new int[num_elems];

          this->elem_orders = new int[num_elems];

          for (unsigned int coeffs_i = 0; coeffs_i < num_coeffs; coeffs_i++)
            this->mono_coeffs[parsed_xml_solution->mono_coeffs().at(coeffs_i).id()] = std::complex<double>(parsed_xml_solution->mono_coeffs().at(coeffs_i).re(), parsed_xml_solution->mono_coeffs().at(coeffs_i).im().get());

          for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
            this->elem_orders[parsed_xml_solution->elem_orders().at(elems_i).id()] = parsed_xml_solution->elem_orders().at(elems_i).ord();

          for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
            for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
              this->elem_coeffs[component_i][parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).id()] = parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).c();
        }

        init_dxdy_buffer();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionLoadFailureException(e.what());
      }
    }

    template<typename Scalar>
    Scalar Solution<Scalar>::get_ref_value(Element* e, double xi1, double xi2, int component, int item)
    {
      if(e==NULL) 
        throw Exceptions::NullException(1);

      set_active_element(e);

      int o = elem_orders[e->id];
      Scalar* mono = dxdy_coeffs[component][item];
      Scalar result = 0.0;
      int k = 0;
      for (int i = 0; i <= o; i++)
      {
        Scalar row = mono[k++];
        for (int j = 0; j < (this->mode ? o : i); j++)
          row = row * xi1 + mono[k++];
        result = result * xi2 + row;
      }
      return result;
    }

    template<typename Scalar>
    Scalar Solution<Scalar>::get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b)
    {
      if(e == NULL) 
        throw Exceptions::NullException(1);

      if(this->num_components == 1)
      {
        if(b == 0)
          return get_ref_value(e, xi1, xi2, a, b);
        if(b == 1 || b == 2)
        {
          double2x2 m;
          double xx, yy;
          this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
          Scalar dx = get_ref_value(e_last = e, xi1, xi2, a, 1);
          Scalar dy = get_ref_value(e, xi1, xi2, a, 2);
          if(b == 1) return m[0][0]*dx + m[0][1]*dy; // H2D_FN_DX
          if(b == 2) return m[1][0]*dx + m[1][1]*dy; // H2D_FN_DY
        }
        else
        {
#ifdef H2D_USE_SECOND_DERIVATIVES
          double2x2 mat;
          double3x2 mat2;
          double xx, yy;

          this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, mat);
          this->refmap->second_ref_map_at_point(xi1, xi2, xx, yy, mat2);

          Scalar vx = get_ref_value(e, xi1, xi2, a, 1);
          Scalar vy = get_ref_value(e, xi1, xi2, a, 2);
          Scalar vxx = get_ref_value(e, xi1, xi2, a, 3);
          Scalar vyy = get_ref_value(e, xi1, xi2, a, 4);
          Scalar vxy = get_ref_value(e, xi1, xi2, a, 5);
          if(b == 3)
            return sqr(mat[0][0])*vxx + 2*mat[0][1]*mat[0][0]*vxy + sqr(mat[0][1])*vyy + mat2[0][0]*vx + mat2[0][1]*vy;   // dxx
          if(b == 4)
            return sqr(mat[1][0])*vxx + 2*mat[1][1]*mat[1][0]*vxy + sqr(mat[1][1])*vyy + mat2[2][0]*vx + mat2[2][1]*vy;   // dyy
          if(b == 5)
            return mat[0][0]*mat[1][0]*vxx + (mat[0][0]*mat[1][1] + mat[1][0]*mat[0][1])*vxy + mat[0][1]*mat[1][1]*vyy + mat2[1][0]*vx + mat2[1][1]*vy;   //dxy
#endif
        }
      }
      else // vector solution
      {
        if(b == 0)
        {
          double2x2 m;
          double xx, yy;
          this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
          Scalar vx = get_ref_value(e, xi1, xi2, 0, 0);
          Scalar vy = get_ref_value(e, xi1, xi2, 1, 0);
          if(a == 0) return m[0][0]*vx + m[0][1]*vy; // H2D_FN_VAL_0
          if(a == 1) return m[1][0]*vx + m[1][1]*vy; // H2D_FN_VAL_1
        }
        else
          throw Hermes::Exceptions::Exception("Getting derivatives of the vector solution: Not implemented yet.");
      }
      throw Hermes::Exceptions::Exception("internal error: reached end of non-void function");
      return 0;
    }

    template<typename Scalar>
    Scalar** Solution<Scalar>::get_ref_values_transformed(Element* e, double x, double y)
    {
      set_active_element(e);

      double x_ref, y_ref;
      double x_dummy, y_dummy;

      this->get_refmap()->untransform(e, x, y, x_ref, y_ref);

      Scalar** toReturn = new Scalar*[2];
      double2x2 mat;
      double3x2 mat2;
      this->refmap->inv_ref_map_at_point(x_ref, y_ref, x_dummy, y_dummy, mat);
      this->refmap->second_ref_map_at_point(x_ref, y_ref, x_dummy, y_dummy, mat2);

      if(this->num_components == 1)
      {
        toReturn[0] = new Scalar[6];

        int o = elem_orders[e->id];

        Scalar result[6];
#ifdef H2D_USE_SECOND_DERIVATIVES
        for(int item = 0; item < 6; item++)
#else
        for(int item = 0; item < 3; item++)
#endif
        {
          Scalar* mono = dxdy_coeffs[0][item];
          Scalar result_local = 0.0;
          int k = 0;
          for (int i = 0; i <= o; i++)
          {
            Scalar row = mono[k++];
            for (int j = 0; j < (this->mode ? o : i); j++)
              row = row * x_ref + mono[k++];
            result[item] = result_local * y_ref + row;
          }
        }

        toReturn[0][0] = result[0];
        toReturn[0][1] = mat[0][0]*result[1] + mat[0][1]*result[2];
        toReturn[0][2] = mat[1][0]*result[1] + mat[1][1]*result[2];
#ifdef H2D_USE_SECOND_DERIVATIVES
        toReturn[0][3] = sqr(mat[0][0])*result[3] + 2*mat[0][1]*mat[0][0]*result[5] + sqr(mat[0][1])*result[4] + mat2[0][0]*result[1] + mat2[0][1]*result[2];
        toReturn[0][4] = sqr(mat[1][0])*result[3] + 2*mat[1][1]*mat[1][0]*result[5] + sqr(mat[1][1])*result[4] + mat2[2][0]*result[1] + mat2[2][1]*result[2];
        toReturn[0][5] = mat[0][0]*mat[1][0]*result[3] + (mat[0][0]*mat[1][1] + mat[1][0]*mat[0][1])*result[5] + mat[0][1]*mat[1][1]*result[4] + mat2[1][0]*result[1] + mat2[1][1]*result[2];
#endif
      }
      else // vector solution
      {
        toReturn[0] = new Scalar[1];

        Scalar vx = get_ref_value(e, x_ref, y_ref, 0, 0);
        Scalar vy = get_ref_value(e, x_ref, y_ref, 1, 0);
        toReturn[0][0] = mat[0][0]*vx + mat[0][1]*vy;
        toReturn[1][0] = mat[1][0]*vx + mat[1][1]*vy;
      }

      return toReturn;
    }

    template<typename Scalar>
    Func<Scalar>* Solution<Scalar>::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
    {
      double xi1, xi2;

      Func<Scalar>* toReturn = new Func<Scalar>(1, this->num_components);

      if(sln_type == HERMES_EXACT)
      {
        if(this->num_components == 1)
        {
          Scalar val, dx = 0.0, dy = 0.0;
          val = (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_function(x, y, dx, dy);
          toReturn->val = new Scalar[1];
          toReturn->dx = new Scalar[1];
          toReturn->dy = new Scalar[1];
          toReturn->val[0] = val * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dx[0] = dx * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dy[0] = dy * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
        }
        else
        {
          Scalar2<Scalar> dx(0.0, 0.0), dy(0.0, 0.0);
          Scalar2<Scalar> val = (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_function(x, y, dx, dy);

          toReturn->val0 = new Scalar[1];
          toReturn->dx0 = new Scalar[1];
          toReturn->dy0 = new Scalar[1];
          toReturn->val1 = new Scalar[1];
          toReturn->dx1 = new Scalar[1];
          toReturn->dy1 = new Scalar[1];
          toReturn->val0[0] = val[0] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->val1[0] = val[1] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dx0[0] = dx[0] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dx1[0] = dx[1] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dy0[0] = dy[0] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          toReturn->dy1[0] = dy[1] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
        }
#ifdef H2D_USE_SECOND_DERIVATIVES
        this->warn("Cannot obtain second derivatives of an exact solution.");
#endif
        return toReturn;
      }
      else if(sln_type == HERMES_UNDEF)
      {
        throw Hermes::Exceptions::Exception("Cannot obtain values -- uninitialized solution. The solution was either "
          "not calculated yet or you used the assignment operator which destroys "
          "the solution on its right-hand side.");
        return NULL;
      }
      else // HERMES_SLN
      {
        if(e == NULL)
          e = RefMap::element_on_physical_coordinates(use_MeshHashGrid, this->mesh, x, y, &xi1, &xi2);
        else
          RefMap::untransform(e, x, y, xi1, xi2);

        if(e != NULL)
        {
          if(this->num_components == 1)
          {
            toReturn->val = new Scalar[1];
            toReturn->dx = new Scalar[1];
            toReturn->dy = new Scalar[1];

            toReturn->val[0] = get_ref_value(e, xi1, xi2, 0, 0);

            double2x2 m;
            double xx, yy;
            this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
            Scalar dx = get_ref_value(e, xi1, xi2, 0, 1);
            Scalar dy = get_ref_value(e, xi1, xi2, 0, 2);
            toReturn->dx[0] = m[0][0]*dx + m[0][1]*dy;
            toReturn->dy[0] = m[1][0]*dx + m[1][1]*dy;

#ifdef H2D_USE_SECOND_DERIVATIVES
            toReturn->laplace = new Scalar[1];
            double2x2 mat;
            double3x2 mat2;

            this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, mat);
            this->refmap->second_ref_map_at_point(xi1, xi2, xx, yy, mat2);

            Scalar vxx = get_ref_value(e, xi1, xi2, 0, 3);
            Scalar vyy = get_ref_value(e, xi1, xi2, 0, 4);
            Scalar vxy = get_ref_value(e, xi1, xi2, 0, 5);
            Scalar dxx = sqr(mat[0][0])*vxx + 2*mat[0][1]*mat[0][0]*vxy + sqr(mat[0][1])*vyy + mat2[0][0]*dx + mat2[0][1]*dy;   // dxx
            Scalar dyy = sqr(mat[1][0])*vxx + 2*mat[1][1]*mat[1][0]*vxy + sqr(mat[1][1])*vyy + mat2[2][0]*dx + mat2[2][1]*dy;   // dyy
            toReturn->laplace[0] = dxx + dyy;
#endif
          }
          else // vector solution
          {
            toReturn->val0 = new Scalar[1];
            toReturn->val1 = new Scalar[1];

            double2x2 m;
            double xx, yy;
            this->refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
            Scalar vx = get_ref_value(e, xi1, xi2, 0, 0);
            Scalar vy = get_ref_value(e, xi1, xi2, 1, 0);
            toReturn->val0[0] = m[0][0]*vx + m[0][1]*vy;
            toReturn->val1[0] = m[1][0]*vx + m[1][1]*vy;
            Hermes::Mixins::Loggable::Static::warn("Derivatives of vector functions not implemented yet.");
          }
          return toReturn;
        }

        this->warn("Point (%g, %g) does not lie in any element.", x, y);
        return NULL;
      }
    }

    template class HERMES_API Solution<double>;
    template class HERMES_API Solution<std::complex<double> >;
  }
}
