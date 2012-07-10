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

#include "solution_h2d_xml.h"

#include <iostream>
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
    void Solution<Scalar>::init()
    {
      memset(tables, 0, sizeof(tables));
      memset(elems,  0, sizeof(elems));
      memset(oldest, 0, sizeof(oldest));
      transform = true;
      sln_type = HERMES_UNDEF;
      sln_vector = NULL;
      space = NULL;
      this->num_components = 0;
      e_last = NULL;

      for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
          tables[i][j] = new std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>;

      mono_coeffs = NULL;
      elem_coeffs[0] = elem_coeffs[1] = NULL;
      elem_orders = NULL;
      dxdy_buffer = NULL;
      num_coeffs = num_elems = 0;
      num_dofs = -1;

      this->set_quad_2d(&g_quad_2d_std);
    }

    template<typename Scalar>
    Solution<Scalar>::Solution()
        : MeshFunction<Scalar>()
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
    }

    template<typename Scalar>
    Solution<Scalar>::Solution(const Mesh *mesh) : MeshFunction<Scalar>(mesh)
    {
      space_type = HERMES_INVALID_SPACE;
      this->init();
      this->mesh = mesh;
    }

    template<typename Scalar>
    Solution<Scalar>::Solution(Space<Scalar>* s, Vector<Scalar>* coeff_vec) : MeshFunction<Scalar>(s->get_mesh())
    {
      space_type = s->get_type();
      space = s;
      this->init();
      this->mesh = s->get_mesh();
      Solution<Scalar>::vector_to_solution(coeff_vec, s, this);
    }

    template<typename Scalar>
    Solution<Scalar>::Solution(Space<Scalar>* s, Scalar* coeff_vec) : MeshFunction<Scalar>(s->get_mesh())
    {
      space_type = s->get_type();
      space = s;
      this->init();
      this->mesh = s->get_mesh();
      Solution<Scalar>::vector_to_solution(coeff_vec, s, this);
    }

    template<typename Scalar>
    void Solution<Scalar>::assign(Solution<Scalar>* sln)
    {
      if(sln->sln_type == HERMES_UNDEF) throw Hermes::Exceptions::Exception("Solution being assigned is uninitialized.");
      if(sln->sln_type != HERMES_SLN) { copy(sln); return; }

      free();

      this->mesh = sln->mesh;
      // Solution vector and space setting.
      this->sln_vector = sln->sln_vector;
      space = sln->space;
      space_type = sln->get_space_type();
      space_seq = sln->get_space_seq();

      mono_coeffs = sln->mono_coeffs;        sln->mono_coeffs = NULL;
      elem_coeffs[0] = sln->elem_coeffs[0];  sln->elem_coeffs[0] = NULL;
      elem_coeffs[1] = sln->elem_coeffs[1];  sln->elem_coeffs[1] = NULL;
      elem_orders = sln->elem_orders;      sln->elem_orders = NULL;
      dxdy_buffer = sln->dxdy_buffer;      sln->dxdy_buffer = NULL;
      num_coeffs = sln->num_coeffs;          sln->num_coeffs = 0;
      num_elems = sln->num_elems;          sln->num_elems = 0;

      sln_type = sln->sln_type;
      this->num_components = sln->num_components;

      memset(sln->tables, 0, sizeof(sln->tables));
    }

    template<typename Scalar>
    void Solution<Scalar>::copy(const Solution<Scalar>* sln)
    {
      if(sln->sln_type == HERMES_UNDEF) throw Hermes::Exceptions::Exception("Solution being copied is uninitialized.");

      free();

      this->mesh = sln->mesh;

      sln_type = sln->sln_type;
      space_type = sln->get_space_type();
      this->num_components = sln->num_components;
      num_dofs = sln->num_dofs;

      if(sln->sln_type == HERMES_SLN) // standard solution: copy coefficient arrays
      {
        num_coeffs = sln->num_coeffs;
        num_elems = sln->num_elems;

        mono_coeffs = new Scalar[num_coeffs];
        memcpy(mono_coeffs, sln->mono_coeffs, sizeof(Scalar) * num_coeffs);

        for (int l = 0; l < this->num_components; l++)
        {
          elem_coeffs[l] = new int[num_elems];
          memcpy(elem_coeffs[l], sln->elem_coeffs[l], sizeof(int) * num_elems);
        }

        elem_orders = new int[num_elems];
        memcpy(elem_orders, sln->elem_orders, sizeof(int) * num_elems);

        init_dxdy_buffer();

        if(this->sln_vector != NULL)
          delete [] this->sln_vector;
        this->sln_vector = new Scalar[sln->space->get_num_dofs()];
        for(int i = 0; i < sln->space->get_num_dofs(); i++)
          this->sln_vector[i] = sln->sln_vector[i];
      }
      else // Const, exact handled differently.
        throw Hermes::Exceptions::Exception("Undefined or exact solutions cannot be copied into an instance of Solution already coming from computation.");

      space = sln->space;
      space_seq = sln->space_seq;

      this->element = NULL;
    }

    template<typename Scalar>
    MeshFunction<Scalar>* Solution<Scalar>::clone()
    {
      Solution<Scalar>* sln = new Solution<Scalar>();
      sln->copy(this);
      return sln;
    }

    template<typename Scalar>
    void Solution<Scalar>::free_tables()
    {
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          if(tables[i][j] != NULL)
          {
            for(typename std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>::iterator it = tables[i][j]->begin(); it != tables[i][j]->end(); it++)
            {
              for(unsigned int l = 0; l < it->second->get_size(); l++)
                if(it->second->present(l))
                  ::free(it->second->get(l));
              delete it->second;
            }
            delete tables[i][j];
            tables[i][j] = NULL;
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

        space = NULL;

        if(this->sln_vector != NULL)
        {
          delete [] this->sln_vector;
          this->sln_vector = NULL;
        }
    }

    template<typename Scalar>
    Solution<Scalar>::~Solution()
    {
      free();
      space_type = HERMES_INVALID_SPACE;
      space = NULL;
    }

    static struct mono_lu_init
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
    void Solution<Scalar>::set_coeff_vector(const Space<Scalar>* space, const Vector<Scalar>* vec,
        bool add_dir_lift, int start_index)
    {
      // Sanity check.
      if(space == NULL) throw Exceptions::NullException(1);
      if(vec == NULL) throw Exceptions::NullException(2);

      space_type = space->get_type();
      Scalar* coeffs = new Scalar[vec->length()];
      vec->extract(coeffs);
      this->set_coeff_vector(space, coeffs, add_dir_lift, start_index);
      delete [] coeffs;
    }

    template<typename Scalar>
    void Solution<Scalar>::set_coeff_vector(const Space<Scalar>* space, const Scalar* coeffs,
        bool add_dir_lift, int start_index)
    {
      // Sanity check.
      if(space == NULL) throw Exceptions::NullException(1);

      // Initialize precalc shapeset using the space's shapeset.
      Shapeset *shapeset = space->shapeset;
      if(space->shapeset == NULL) throw Exceptions::Exception("Space->shapeset == NULL in Solution<Scalar>::set_coeff_vector().");
      PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
      if(pss == NULL) throw Exceptions::Exception("PrecalcShapeset could not be allocated in Solution<Scalar>::set_coeff_vector().");
      set_coeff_vector(space, pss, coeffs, add_dir_lift, start_index);
      delete pss;
    }

    template<typename Scalar>
    void Solution<Scalar>::set_coeff_vector(const Space<Scalar>* space, PrecalcShapeset* pss,
        const Scalar* coeff_vec, bool add_dir_lift, int start_index)
    {
      int o;

      // Sanity checks.
      if(space == NULL) throw Exceptions::NullException(1);
      if(space->get_mesh() == NULL) throw Exceptions::Exception("Mesh == NULL in Solution<Scalar>::set_coeff_vector().");
      if(pss == NULL) throw Exceptions::NullException(2);
      if(coeff_vec == NULL) throw Exceptions::NullException(3);
      if(coeff_vec == NULL) throw Exceptions::Exception("Coefficient vector == NULL in Solution<Scalar>::set_coeff_vector().");
      if(!space->is_up_to_date())
        throw Exceptions::Exception("Provided 'space' is not up to date.");
      if(space->shapeset != pss->shapeset)
        throw Exceptions::Exception("Provided 'space' and 'pss' must have the same shapesets.");

      free();

      if(this->sln_vector != NULL)
        delete [] this->sln_vector;
      int ndof = space->get_num_dofs();
      this->sln_vector = new Scalar[ndof];
      for(int i = 0; i < ndof; i++)
      {
        // By adding start_index we move to the desired section of coeff_vec.
        this->sln_vector[i] = coeff_vec[i + start_index];
      }

      this->space_type = space->get_type();
      this->space = space;
      this->space_seq = space->get_seq();

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
        for (unsigned int k = 0; k < e->get_num_surf(); k++)
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

      if(this->mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL.\n");
      init_dxdy_buffer();
      this->element = NULL;
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Scalar* solution_vector,
        Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> solutions,
        Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector == NULL) throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

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

      for(unsigned int i = 0; i < solutions.size(); i++)
      {
        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }

      return;
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Scalar* solution_vector, const Space<Scalar>* space,
        Solution<Scalar>* solution, bool add_dir_lift, int start_index)
    {
      // Sanity checks.
      if(solution_vector == NULL) throw Exceptions::NullException(1);
      if(space == NULL) throw Exceptions::NullException(2);
      if(solution == NULL) throw Exceptions::NullException(3);

      solution->set_coeff_vector(space, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Vector<Scalar>* solution_vector, Hermes::vector<const Space<Scalar>*> spaces,
      Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector == NULL) throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

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

      for(unsigned int i = 0; i < solutions.size(); i++)
      {
        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }

      return;
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Vector<Scalar>* solution_vector, const Space<Scalar>* space,
        Solution<Scalar>* solution, bool add_dir_lift, int start_index)
    {
      // Sanity checks.
        if(solution_vector == NULL) throw Exceptions::NullException(1);
      if(space == NULL) throw Exceptions::NullException(2);
      if(solution == NULL) throw Exceptions::NullException(3);

      solution->set_coeff_vector(space, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solutions(const Scalar* solution_vector, Hermes::vector<const Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<PrecalcShapeset *> pss,
        Hermes::vector<bool> add_dir_lift, Hermes::vector<int> start_indices)
    {
      if(solution_vector==NULL) throw Exceptions::NullException(1);
      if(spaces.size() != solutions.size()) throw Exceptions::LengthException(2, 3, spaces.size(), solutions.size());

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

      for(unsigned int i = 0; i < solutions.size(); i++)
      {
        if(add_dir_lift == Hermes::vector<bool>())
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], true, start_indices_new[i]);
        else
          Solution<Scalar>::vector_to_solution(solution_vector, spaces[i], solutions[i], add_dir_lift[i], start_indices_new[i]);
      }

      return;
    }

    template<typename Scalar>
    void Solution<Scalar>::vector_to_solution(const Scalar* solution_vector, const Space<Scalar>* space, Solution<Scalar>* solution,
        PrecalcShapeset* pss, bool add_dir_lift, int start_index)
    {
      if(solution_vector == NULL) throw Exceptions::NullException(1);
      if(space == NULL) throw Exceptions::NullException(2);
      if(solution == NULL) throw Exceptions::NullException(3);
      if(pss == NULL) throw Exceptions::NullException(4);

      solution->set_coeff_vector(space, pss, solution_vector, add_dir_lift, start_index);
    }

    template<typename Scalar>
    void Solution<Scalar>::set_dirichlet_lift(const Space<Scalar>* space, PrecalcShapeset* pss)
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
      if(transform != enable) free_tables();
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
        dynamic_cast<ExactSolution<Scalar>*>(this)->exact_multiplicator *= coef;
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
      // if(e == element) return; // FIXME
      if(!e->active) throw Hermes::Exceptions::Exception("Cannot select inactive element. Wrong mesh?");
      MeshFunction<Scalar>::set_active_element(e);

      // try finding an existing table for e
      for (cur_elem = 0; cur_elem < 4; cur_elem++)
        if(elems[this->cur_quad][cur_elem] == e)
          break;

      // if not found, free the oldest one and use its slot
      if(cur_elem >= 4)
      {
        if(tables[this->cur_quad][oldest[this->cur_quad]] != NULL)
        {
          for(typename std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>::iterator it = tables[this->cur_quad][oldest[this->cur_quad]]->begin(); it != tables[this->cur_quad][oldest[this->cur_quad]]->end(); it++)
          {
            for(unsigned int l = 0; l < it->second->get_size(); l++)
              if(it->second->present(l))
                ::free(it->second->get(l));
            delete it->second;
          }
          delete tables[this->cur_quad][oldest[this->cur_quad]];
          tables[this->cur_quad][oldest[this->cur_quad]] = NULL;
          elems[this->cur_quad][oldest[this->cur_quad]] = NULL;
        }

        tables[this->cur_quad][oldest[this->cur_quad]] = new std::map<uint64_t, LightArray<struct Function<Scalar>::Node*>*>;

        cur_elem = oldest[this->cur_quad];
        if(++oldest[this->cur_quad] >= 4)
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

      this->sub_tables = tables[this->cur_quad][cur_elem];

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
      double3x2 *mat2, *mm;
      int i, mstep = 0;

      // H1 space
      if(space_type == HERMES_H1_SPACE)
      {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
        if(((newmask & H2D_SECOND) == H2D_SECOND && (oldmask & H2D_SECOND) != H2D_SECOND))
        {
          this->update_refmap();
          mat = this->refmap->get_inv_ref_map(order);
          mat2 = this->refmap->get_second_ref_map(order);
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
      if(sln_type == HERMES_EXACT)
        throw Exceptions::Exception("Exact solution cannot be saved to a file.");
      if(sln_type == HERMES_UNDEF)
        throw Exceptions::Exception("Cannot save -- uninitialized solution.");

      try
      {
        XMLSolution::solution xmlsolution(XMLSolution::sln_vector(), this->num_components,
          this->num_elems, this->num_coeffs, this->space->get_num_dofs());

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

        for(unsigned int sln_coeff_i = 0; sln_coeff_i < this->space->get_num_dofs(); sln_coeff_i++)
          xmlsolution.sln_vector().sln_coeff().push_back(XMLSolution::sln_coeff(sln_coeff_i, this->sln_vector[sln_coeff_i]));

        std::string solution_schema_location(H2D_XML_SCHEMAS_DIRECTORY);
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
      return;
    }

    template<>
    void Solution<std::complex<double> >::save(const char* filename) const
    {
      if(sln_type == HERMES_EXACT)
        throw Exceptions::Exception("Exact solution cannot be saved to a file.");
      if(sln_type == HERMES_UNDEF)
        throw Exceptions::Exception("Cannot save -- uninitialized solution.");

      try
      {
        XMLSolution::solution xmlsolution(XMLSolution::sln_vector(), this->num_components, this->num_elems, this->num_coeffs, this->space->get_num_dofs());

        for(unsigned int coeffs_i = 0; coeffs_i < this->num_coeffs; coeffs_i++)
        {
          xmlsolution.mono_coeffs().push_back(XMLSolution::mono_coeffs(coeffs_i, mono_coeffs[coeffs_i].real()));
          xmlsolution.mono_coeffs().back().imaginary() = mono_coeffs[coeffs_i].imag();
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

        for(unsigned int sln_coeff_i = 0; sln_coeff_i < this->space->get_num_dofs(); sln_coeff_i++)
        {
          xmlsolution.sln_vector().sln_coeff().push_back(XMLSolution::sln_coeff(sln_coeff_i, this->sln_vector[sln_coeff_i].real()));
          xmlsolution.sln_vector().sln_coeff().back().imaginary() = this->sln_vector[sln_coeff_i].imag();
        }

        std::string solution_schema_location(H2D_XML_SCHEMAS_DIRECTORY);
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
      return;
    }

    template<>
    void Solution<double>::load(const char* filename, Mesh* mesh)
    {
      free();
      sln_type = HERMES_SLN;
      this->mesh = mesh;

      try
      {
        std::auto_ptr<XMLSolution::solution> parsed_xml_solution(XMLSolution::solution_(filename));

        this->num_coeffs = parsed_xml_solution->num_coeffs();
        this->num_elems = parsed_xml_solution->num_elems();
        this->num_components = parsed_xml_solution->num_components();
        this->num_dofs = parsed_xml_solution->num_dofs();

        this->mono_coeffs = new double[num_coeffs];
        memset(this->mono_coeffs, 0, this->num_coeffs*sizeof(double));

        for(unsigned int component_i = 0; component_i < num_components; component_i++)
          elem_coeffs[component_i] = new int[num_elems];

        this->elem_orders = new int[num_elems];
        this->sln_vector = new double[num_dofs];

        for (unsigned int coeffs_i = 0; coeffs_i < num_coeffs; coeffs_i++)
          this->mono_coeffs[parsed_xml_solution->mono_coeffs().at(coeffs_i).id()] = parsed_xml_solution->mono_coeffs().at(coeffs_i).real();

        for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
          this->elem_orders[parsed_xml_solution->elem_orders().at(elems_i).id()] = parsed_xml_solution->elem_orders().at(elems_i).order();

        for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
          for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
            this->elem_coeffs[component_i][parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).id()] = parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).coeff();

        for(unsigned int sln_coeff_i = 0; sln_coeff_i < this->num_dofs; sln_coeff_i++)
          this->sln_vector[parsed_xml_solution->sln_vector().sln_coeff().at(sln_coeff_i).id()] = parsed_xml_solution->sln_vector().sln_coeff().at(sln_coeff_i).real();

        init_dxdy_buffer();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionLoadFailureException(e.what());
      }
      return;
    }

    template<>
    void Solution<std::complex<double> >::load(const char* filename, Mesh* mesh)
    {
      free();
      sln_type = HERMES_SLN;
      this->mesh = mesh;

      try
      {
        std::auto_ptr<XMLSolution::solution> parsed_xml_solution(XMLSolution::solution_(filename));

        this->num_coeffs = parsed_xml_solution->num_coeffs();
        this->num_elems = parsed_xml_solution->num_elems();
        this->num_components = parsed_xml_solution->num_components();
        this->num_dofs = parsed_xml_solution->num_dofs();

        this->mono_coeffs = new std::complex<double>[num_coeffs];
        memset(this->mono_coeffs, 0, this->num_coeffs*sizeof(std::complex<double>));

        for(unsigned int component_i = 0; component_i < num_components; component_i++)
          elem_coeffs[component_i] = new int[num_elems];

        this->elem_orders = new int[num_elems];
        this->sln_vector = new std::complex<double>[num_dofs];

        for (unsigned int coeffs_i = 0; coeffs_i < num_coeffs; coeffs_i++)
          this->mono_coeffs[parsed_xml_solution->mono_coeffs().at(coeffs_i).id()] = std::complex<double>(parsed_xml_solution->mono_coeffs().at(coeffs_i).real(), parsed_xml_solution->mono_coeffs().at(coeffs_i).imaginary().get());

        for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
          this->elem_orders[parsed_xml_solution->elem_orders().at(elems_i).id()] = parsed_xml_solution->elem_orders().at(elems_i).order();

        for (unsigned int component_i = 0; component_i < this->num_components; component_i++)
          for (unsigned int elems_i = 0; elems_i < num_elems; elems_i++)
            this->elem_coeffs[component_i][parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).id()] = parsed_xml_solution->component().at(component_i).elem_coeffs().at(elems_i).coeff();

        for(unsigned int sln_coeff_i = 0; sln_coeff_i < this->num_dofs; sln_coeff_i++)
          this->sln_vector[parsed_xml_solution->sln_vector().sln_coeff().at(sln_coeff_i).id()] = std::complex<double>(parsed_xml_solution->sln_vector().sln_coeff().at(sln_coeff_i).real(), parsed_xml_solution->sln_vector().sln_coeff().at(sln_coeff_i).imaginary().get());

        init_dxdy_buffer();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionLoadFailureException(e.what());
      }
      return;
    }

    template<typename Scalar>
    Scalar Solution<Scalar>::get_ref_value(Element* e, double xi1, double xi2, int component, int item)
    {
      if(e==NULL) throw Exceptions::NullException(1);
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

    static bool is_in_ref_domain(Element* e, double xi1, double xi2)
    {
      const double TOL = 1e-11;
      if(e->get_num_surf() == 3)
        return (xi1 + xi2 <= TOL) && (xi1 + 1.0 >= -TOL) && (xi2 + 1.0 >= -TOL);
      else
        return (xi1 - 1.0 <= TOL) && (xi1 + 1.0 >= -TOL) && (xi2 - 1.0 <= TOL) && (xi2 + 1.0 >= -TOL);
    }

    template<typename Scalar>
    Scalar Solution<Scalar>::get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b)
    {
      if(e==NULL) throw Exceptions::NullException(1);
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
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
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
    Scalar Solution<Scalar>::get_pt_value(double x, double y, int item)
    {
      double xi1, xi2;

      int a = 0, b = 0, mask = item; // a = component, b = val, dx, dy, dxx, dyy, dxy
      if(this->num_components == 1) mask = mask & H2D_FN_COMPONENT_0;
      if((mask & (mask - 1)) != 0) throw Hermes::Exceptions::Exception("'item' is invalid. ");
      if(mask >= 0x40) { a = 1; mask >>= 6; }
      while (!(mask & 1)) { mask >>= 1; b++; }

      if(sln_type == HERMES_EXACT)
      {
        if(this->num_components == 1)
        {
          Scalar val, dx = 0.0, dy = 0.0;
          val = (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_function(x, y, dx, dy);
          if(b == 0) return val * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          if(b == 1) return dx * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          if(b == 2) return dy * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
        }
        else
        {
          Scalar2<Scalar> dx(0.0, 0.0), dy(0.0, 0.0);
          Scalar2<Scalar> val = (static_cast<ExactSolutionVector<Scalar>*>(this))->exact_function(x, y, dx, dy);
          if(b == 0) return val[a] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          if(b == 1) return dx[a] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
          if(b == 2) return dy[a] * (static_cast<ExactSolutionScalar<Scalar>*>(this))->exact_multiplicator;
        }
        throw Hermes::Exceptions::Exception("Cannot obtain second derivatives of an exact solution.");
      }
      else if(sln_type == HERMES_UNDEF)
      {
        throw Hermes::Exceptions::Exception("Cannot obtain values -- uninitialized solution. The solution was either "
          "not calculated yet or you used the assignment operator which destroys "
          "the solution on its right-hand side.");
      }

      // try the last visited element and its neighbours
      if(e_last != NULL)
      {
        Element* elem[5];
        elem[0] = e_last;
        for (unsigned int i = 1; i <= e_last->get_num_surf(); i++)
          elem[i] = e_last->get_neighbor(i-1);

        for (unsigned int i = 0; i <= e_last->get_num_surf(); i++)
          if(elem[i] != NULL)
          {
            this->refmap->set_active_element(elem[i]);
            this->refmap->untransform(elem[i], x, y, xi1, xi2);
            if(is_in_ref_domain(elem[i], xi1, xi2))
            {
              e_last = elem[i];
              return get_ref_value_transformed(elem[i], xi1, xi2, a, b);
            }
          }
      }

      // go through all elements
      Element *e;
      for_all_active_elements(e, this->mesh)
      {
        this->refmap->set_active_element(e);
        this->refmap->untransform(e, x, y, xi1, xi2);
        if(is_in_ref_domain(e, xi1, xi2))
        {
          e_last = e;
          return get_ref_value_transformed(e, xi1, xi2, a, b);
        }
      }

      this->warn("Point (%g, %g) does not lie in any element.", x, y);
      return NAN;
    }

    template<typename Scalar>
    const Space<Scalar>* Solution<Scalar>::get_space()
    {
      if(this->sln_type == HERMES_SLN)
        return space;
      else
      {
        throw Hermes::Exceptions::Exception("Solution<Scalar>::get_space() called with an instance where FEM space is not defined.");
        return NULL;
      }
    }

    template<typename Scalar>
    int Solution<Scalar>::get_space_seq()
    {
      if(this->sln_type == HERMES_SLN)
        return space_seq;
      else
      {
        throw Hermes::Exceptions::Exception("Solution<Scalar>::get_space_seq() called with an instance where FEM space is not defined.");
        return NULL;
      }
    }

    template<typename Scalar>
    Scalar* Solution<Scalar>::get_sln_vector()
    {
      if(this->sln_type == HERMES_SLN)
        return sln_vector;
      else
      {
        throw Hermes::Exceptions::Exception("Solution<Scalar>::get_sln_vector() called with an instance where FEM space is not defined.");
        return NULL;
      }
    }

    template class HERMES_API Solution<double>;
    template class HERMES_API Solution<std::complex<double> >;
  }
}