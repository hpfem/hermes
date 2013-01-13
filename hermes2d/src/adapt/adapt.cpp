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

#include "umfpack.h"
#include "adapt.h"
#include "hermes2d.h"
#include "global.h"
#include "limit_order.h"
#include "solution.h"
#include "discrete_problem.h"
#include "refmap.h"
#include "quad_all.h"
#include "traverse.h"
#include "refinement_selectors/optimum_selector.h"
#include "matrix.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Adapt<Scalar>::Adapt(Hermes::vector<Space<Scalar>*> spaces,
      Hermes::vector<ProjNormType> proj_norms) :
    spaces(spaces),
      num_act_elems(-1),
      have_errors(false),
      have_coarse_solutions(false),
      have_reference_solutions(false)
    {
      // sanity check
      if(proj_norms.size() > 0 && spaces.size() != proj_norms.size())
        throw Exceptions::LengthException(1, 2, spaces.size(), proj_norms.size());

      this->num = spaces.size();

      // sanity checks
      if((this->num <= 0) || (this->num > H2D_MAX_COMPONENTS)) throw Exceptions::ValueException("components", this->num, 0, H2D_MAX_COMPONENTS);

      // reset values
      memset(errors, 0, sizeof(errors));
      memset(sln, 0, sizeof(sln));
      memset(rsln, 0, sizeof(rsln));
      own_forms = new bool*[H2D_MAX_COMPONENTS];
      for(int i = 0; i < H2D_MAX_COMPONENTS; i++)
      {
        own_forms[i] = new bool[H2D_MAX_COMPONENTS];
        memset(own_forms[i], 0, H2D_MAX_COMPONENTS * sizeof(bool));
      }

      // if norms were not set by the user, set them to defaults
      // according to spaces
      if(proj_norms.size() == 0)
      {
        for (int i = 0; i < this->num; i++)
        {
          switch (spaces[i]->get_type())
          {
          case HERMES_H1_SPACE: proj_norms.push_back(HERMES_H1_NORM); break;
          case HERMES_HCURL_SPACE: proj_norms.push_back(HERMES_HCURL_NORM); break;
          case HERMES_HDIV_SPACE: proj_norms.push_back(HERMES_HDIV_NORM); break;
          case HERMES_L2_SPACE: proj_norms.push_back(HERMES_L2_NORM); break;
          default: throw Hermes::Exceptions::Exception("Unknown space type in Adapt<Scalar>::Adapt().");
          }
        }
      }

      // assign norm weak forms  according to norms selection
      for (int i = 0; i < this->num; i++)
        for (int j = 0; j < this->num; j++)
        {
          error_form[i][j] = NULL;
          norm_form[i][j] = NULL;
        }

        for (int i = 0; i < this->num; i++)
        {
          error_form[i][i] = new MatrixFormVolError(i, i, proj_norms[i]);
          norm_form[i][i] = error_form[i][i];
          own_forms[i][i] = true;
        }
    }

    template<typename Scalar>
    Adapt<Scalar>::Adapt(Space<Scalar>* space, ProjNormType proj_norm) :
    spaces(Hermes::vector<Space<Scalar>*>()),
      num_act_elems(-1),
      have_errors(false),
      have_coarse_solutions(false),
      have_reference_solutions(false)
    {
      if(space == NULL) throw Exceptions::NullException(1);
      spaces.push_back(space);

      this->num = 1;

      // reset values
      memset(errors, 0, sizeof(errors));
      memset(sln, 0, sizeof(sln));
      memset(rsln, 0, sizeof(rsln));
      own_forms = new bool*[H2D_MAX_COMPONENTS];
      for(int i = 0; i < H2D_MAX_COMPONENTS; i++)
      {
        own_forms[i] = new bool[H2D_MAX_COMPONENTS];
        memset(own_forms[i], 0, H2D_MAX_COMPONENTS * sizeof(bool));
      }

      // if norms were not set by the user, set them to defaults
      // according to spaces
      if(proj_norm == HERMES_UNSET_NORM)
      {
        switch (space->get_type())
        {
        case HERMES_H1_SPACE: proj_norm = HERMES_H1_NORM; break;
        case HERMES_HCURL_SPACE: proj_norm = HERMES_HCURL_NORM; break;
        case HERMES_HDIV_SPACE: proj_norm = HERMES_HDIV_NORM; break;
        case HERMES_L2_SPACE: proj_norm = HERMES_L2_NORM; break;
        default: throw Hermes::Exceptions::Exception("Unknown space type in Adapt<Scalar>::Adapt().");
        }
      }

      // assign norm weak forms  according to norms selection
      error_form[0][0] = new MatrixFormVolError(0, 0, proj_norm);
      norm_form[0][0] = error_form[0][0];
      own_forms[0][0] = true;
    }

    template<typename Scalar>
    Adapt<Scalar>::~Adapt()
    {
      for (int i = 0; i < this->num; i++)
        delete [] errors[i];

      // free error_form
      for (int i = 0; i < this->num; i++)
        for (int j = 0; j < this->num; j++)
          if(error_form[i][j] != NULL && own_forms[i][j])
          {
            delete error_form[i][j];
            own_forms[i][j] = false;
          }
    }

    template<typename Scalar>
    bool Adapt<Scalar>::adapt(Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors, double thr, int strat,
      int regularize, double to_be_processed)
    {
      this->tick();
      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      if(!have_errors)
        throw Exceptions::Exception("element errors have to be calculated first, call Adapt<Scalar>::calc_err_est().");

      if(refinement_selectors.empty())
        throw Exceptions::NullException(1);
      if(spaces.size() != refinement_selectors.size())
        throw Exceptions::LengthException(1, refinement_selectors.size(), spaces.size());

      //get meshes
      int max_id = -1;
      Mesh* meshes[H2D_MAX_COMPONENTS];
      for (int j = 0; j < this->num; j++)
      {
        meshes[j] = this->spaces[j]->get_mesh();
        if(rsln[j] != NULL)
        {
          rsln[j]->set_quad_2d(&g_quad_2d_std);
          rsln[j]->enable_transform(false);
        }
        if(meshes[j]->get_max_element_id() > max_id)
          max_id = meshes[j]->get_max_element_id();
      }

      //reset element refinement info
      int** idx = new int*[max_id];
      for(int i = 0; i < max_id; i++)
        idx[i] = new int[num];

      Element* e;
      for(int j = 0; j < max_id; j++)
        for(int l = 0; l < this->num; l++)
          idx[j][l] = -1; // element not refined

      double err0_squared = 1000.0;
      double processed_error_squared = 0.0;

      std::vector<ElementToRefine> elem_inx_to_proc; //list of indices of elements that are going to be processed
      elem_inx_to_proc.reserve(num_act_elems);

      //adaptivity loop
      double error_squared_threshold = -1; //an error threshold that breaks the adaptivity loop in a case of strategy 1
      int num_ignored_elem = 0; //a number of ignored elements
      int num_not_changed = 0; //a number of element that were not changed
      int num_priority_elem = 0; //a number of elements that were processed using priority queue

      // Structures traversed in reality using strategies.
      Hermes::vector<int> ids;
      Hermes::vector<int> components;
      Hermes::vector<int> current_orders;
      bool first_regular_element = true; // true if first regular element was not processed yet
      bool error_level_reached = false;

      for(int inx_regular_element = 0; inx_regular_element < num_act_elems || !priority_queue.empty();)
      {
        int id, comp;

        // Process the queuse(s) to see what elements to really refine.
        if(priority_queue.empty())
        {
          id = regular_queue[inx_regular_element].id;
          comp = regular_queue[inx_regular_element].comp;
          inx_regular_element++;

          // Get info linked with the element
          double err_squared = errors[comp][id];

          if(first_regular_element)
          {
            error_squared_threshold = thr * err_squared;
            first_regular_element = false;
          }

          // first refinement strategy:
          // refine elements until prescribed amount of error is processed
          // if more elements have similar error refine all to keep the mesh symmetric
          if((strat == 0) && (processed_error_squared > sqrt(thr) * errors_squared_sum)
            && fabs((err_squared - err0_squared)/err0_squared) > 1e-3)
            error_level_reached = true;

          // second refinement strategy:
          // refine all elements whose error is bigger than some portion of maximal error
          if((strat == 1) && (err_squared < error_squared_threshold))
            error_level_reached = true;

          if((strat == 2) && (err_squared < thr))
            error_level_reached = true;

          if((strat == 3) && ((err_squared < error_squared_threshold) || (processed_error_squared > 1.5 * to_be_processed)))
            error_level_reached = true;

          // Insert the element only if it complies with the strategy.
          if(!error_level_reached)
          {
            err0_squared = err_squared;
            processed_error_squared += err_squared;
            ids.push_back(id);
            components.push_back(comp);
            current_orders.push_back(this->spaces[comp]->get_element_order(id));
            spaces[comp]->edata[id].changed_in_last_adaptation = true;
          }
          else
            if(priority_queue.empty())
              break;
        }
        // Priority - refine no matter what.
        else
        {
          id = priority_queue.front().id;
          comp = priority_queue.front().comp;
          priority_queue.pop();

          // Insert into appropripate arrays.
          ids.push_back(id);
          components.push_back(comp);
          current_orders.push_back(this->spaces[comp]->get_element_order(id));
          spaces[comp]->edata[id].changed_in_last_adaptation = true;
        }
      }

      if(ids.empty())
      {
        this->warn("None of the elements selected for refinement was refined. Adaptivity step successful, returning 'true'.");
        return true;
      }

      // RefinementSelectors cloning.
      RefinementSelectors::Selector<Scalar>*** global_refinement_selectors = new RefinementSelectors::Selector<Scalar>**[Hermes::Hermes2D::Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];

      for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
      {
        global_refinement_selectors[i] = new RefinementSelectors::Selector<Scalar>*[refinement_selectors.size()];
        for (unsigned int j = 0; j < refinement_selectors.size(); j++)
        {
          if(i == 0)
            global_refinement_selectors[i][j] = refinement_selectors[j];
          else
          {
            global_refinement_selectors[i][j] = refinement_selectors[j]->clone();
            if(dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j]) != NULL)
            {
              dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j])->cached_shape_vals_valid = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(refinement_selectors[j])->cached_shape_vals_valid;
              dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j])->cached_shape_ortho_vals = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(refinement_selectors[j])->cached_shape_ortho_vals;
              dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(global_refinement_selectors[i][j])->cached_shape_vals = dynamic_cast<RefinementSelectors::ProjBasedSelector<Scalar>*>(refinement_selectors[j])->cached_shape_vals;
            }
            if(dynamic_cast<RefinementSelectors::OptimumSelector<Scalar>*>(global_refinement_selectors[i][j]) != NULL)
              dynamic_cast<RefinementSelectors::OptimumSelector<Scalar>*>(global_refinement_selectors[i][j])->num_shapes = dynamic_cast<RefinementSelectors::OptimumSelector<Scalar>*>(refinement_selectors[j])->num_shapes;
          }
        }
      }

      // Solution cloning.
      Solution<Scalar>*** rslns = new Solution<Scalar>**[Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];

      for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
      {
        rslns[i] = new Solution<Scalar>*[this->num];
        for (int j = 0; j < this->num; j++)
        {
          if(rsln[j] != NULL)
            rslns[i][j] = dynamic_cast<Solution<Scalar>*>(rsln[j]->clone());
        }
      }

      this->tick();
      this->info("Adaptivity: data preparation duration: %f s.", this->last());

      // For statistics.
      Hermes::vector<int> numberOfCandidates;


      // The loop
      RefinementSelectors::Selector<Scalar>** current_refinement_selectors;
      Solution<Scalar>** current_rslns;
      int id_to_refine;
#define CHUNKSIZE 1
      int num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);
#pragma omp parallel shared(ids, components, elem_inx_to_proc, meshes, current_orders) private(current_refinement_selectors, current_rslns, id_to_refine) num_threads(num_threads_used)
      {
#pragma omp for schedule(dynamic, CHUNKSIZE)
        for(id_to_refine = 0; id_to_refine < ids.size(); id_to_refine++)
        {
          try
          {
            current_refinement_selectors = global_refinement_selectors[omp_get_thread_num()];
            current_rslns = rslns[omp_get_thread_num()];

            // Get refinement suggestion
            ElementToRefine elem_ref(ids[id_to_refine], components[id_to_refine]);

            // rsln[comp] may be unset if refinement_selectors[comp] == HOnlySelector or POnlySelector
            bool refined = current_refinement_selectors[components[id_to_refine]]->select_refinement(meshes[components[id_to_refine]]->get_element(ids[id_to_refine]), current_orders[id_to_refine], current_rslns[components[id_to_refine]], elem_ref);
            
#pragma omp critical (number_of_candidates)
						{
							if(dynamic_cast<Hermes::Hermes2D::RefinementSelectors::OptimumSelector<Scalar>*>(current_refinement_selectors[components[id_to_refine]]) != NULL)
								numberOfCandidates.push_back(dynamic_cast<Hermes::Hermes2D::RefinementSelectors::OptimumSelector<Scalar>*>(current_refinement_selectors[components[id_to_refine]])->get_candidates().size());
						}

            //add to a list of elements that are going to be refined
  #pragma omp critical (elem_inx_to_proc)
            {
              idx[ids[id_to_refine]][components[id_to_refine]] = (int)elem_inx_to_proc.size();
              elem_inx_to_proc.push_back(elem_ref);
            }
          }
          catch(Hermes::Exceptions::Exception& exception)
          {
            if(this->caughtException == NULL)
              this->caughtException = exception.clone();
          }
          catch(std::exception& exception)
          {
            if(this->caughtException == NULL)
              this->caughtException = new Hermes::Exceptions::Exception(exception.what());
          }
        }
      }

      if(this->caughtException == NULL)
      {
        int averageNumberOfCandidates = 0;
        for(int i = 0; i < numberOfCandidates.size(); i++)
          averageNumberOfCandidates += numberOfCandidates[i];
        averageNumberOfCandidates = averageNumberOfCandidates / numberOfCandidates.size();

        this->info("Adaptivity: total number of refined Elements: %i.", ids.size());
        this->info("Adaptivity: average number of candidates per refined Element: %i.", averageNumberOfCandidates);
      }

      this->tick();
      this->info("Adaptivity: refinement selection duration: %f s.", this->last());

      if(this->caughtException == NULL)
        fix_shared_mesh_refinements(meshes, elem_inx_to_proc, idx, global_refinement_selectors);

      for(unsigned int i = 0; i < Hermes::Hermes2D::Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
      {
        if(i > 0)
          for (unsigned int j = 0; j < refinement_selectors.size(); j++)
            delete global_refinement_selectors[i][j];
        delete [] global_refinement_selectors[i];
      }
      delete [] global_refinement_selectors;

      for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
      {
        if(rslns[i] != NULL)
        {
          for (unsigned int j = 0; j < this->num; j++)
            if(rsln[j] != NULL)
              delete rslns[i][j];
          delete [] rslns[i];
        }
      }
      delete [] rslns;

      for(int i = 0; i < max_id; i++)
        delete [] idx[i];
      delete [] idx;

      if(this->caughtException != NULL)
        throw *(this->caughtException);
      
      //apply refinements
      apply_refinements(elem_inx_to_proc);

      // in singlemesh case, impose same orders across meshes
      homogenize_shared_mesh_orders(meshes);

      // mesh regularization
      if(regularize >= 0)
      {
        if(regularize == 0)
        {
          regularize = 1;
          this->warn("Total mesh regularization is not supported in adaptivity. 1-irregular mesh is used instead.");
        }
        for (int i = 0; i < this->num; i++)
        {
          int* parents;
          parents = meshes[i]->regularize(regularize);
          this->spaces[i]->distribute_orders(meshes[i], parents);
          ::free(parents);
        }
      }

      for (int j = 0; j < this->num; j++)
        if(rsln[j] != NULL)
          rsln[j]->enable_transform(true);

      //store for the user to retrieve
      last_refinements.swap(elem_inx_to_proc);

      have_errors = false;
      if(strat == 2)
        have_errors = true; // space without changes

      // since space changed, assign dofs:
      for(unsigned int i = 0; i < this->spaces.size(); i++)
        this->spaces[i]->assign_dofs();

      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, this->spaces[i]->get_mesh())
          this->spaces[i]->edata[e->id].changed_in_last_adaptation = false;
        for(id_to_refine = 0; id_to_refine < ids.size(); id_to_refine++)
          this->spaces[i]->edata[ids[id_to_refine]].changed_in_last_adaptation = false;
      }

      return false;
    }

    template<typename Scalar>
    Adapt<Scalar>::MatrixFormVolError::MatrixFormVolError(int i, int j) : MatrixFormVol<Scalar>(i, j)
    {
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>* Adapt<Scalar>::MatrixFormVolError::clone() const
    {
      return new MatrixFormVolError(*this);
    }

    template<typename Scalar>
    Adapt<Scalar>::MatrixFormVolError::MatrixFormVolError(int i, int j, ProjNormType type) : MatrixFormVol<Scalar>(i, j), projNormType(type)
    {
    }

    template<typename Scalar>
    template<typename TestFunctionDomain, typename SolFunctionDomain>
    SolFunctionDomain Adapt<Scalar>::MatrixFormVolError::l2_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
      Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
    {
      SolFunctionDomain result = SolFunctionDomain(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]));
      return result;
    }

    template<typename Scalar>
    template<typename TestFunctionDomain, typename SolFunctionDomain>
    SolFunctionDomain Adapt<Scalar>::MatrixFormVolError::h1_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
      Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
    {
      SolFunctionDomain result = SolFunctionDomain(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i])
        + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename Scalar>
    template<typename TestFunctionDomain, typename SolFunctionDomain>
    SolFunctionDomain Adapt<Scalar>::MatrixFormVolError::h1_error_semi_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
      Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
    {
      SolFunctionDomain result = SolFunctionDomain(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename Scalar>
    template<typename TestFunctionDomain, typename SolFunctionDomain>
    SolFunctionDomain Adapt<Scalar>::MatrixFormVolError::hdiv_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
      Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
    {
      throw Hermes::Exceptions::Exception("hdiv error form not implemented yet in hdiv.h.");

      // this is Hcurl code:
      SolFunctionDomain result = SolFunctionDomain(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) +
        u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename Scalar>
    template<typename TestFunctionDomain, typename SolFunctionDomain>
    SolFunctionDomain Adapt<Scalar>::MatrixFormVolError::hcurl_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
      Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
    {
      SolFunctionDomain result = SolFunctionDomain(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) +
        u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename Scalar>
    Scalar Adapt<Scalar>::MatrixFormVolError::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_error_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_NORM:
        return h1_error_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_error_semi_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_error_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_error_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
      default:
        throw Hermes::Exceptions::Exception("Unknown projection type");
        return 0.0;
      }
    }

    template<typename Scalar>
    Hermes::Ord Adapt<Scalar>::MatrixFormVolError::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_error_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_NORM:
        return h1_error_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_error_semi_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_error_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_error_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
      default:
        throw Hermes::Exceptions::Exception("Unknown projection type");
        return Hermes::Ord();
      }
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_est(Solution<Scalar>*sln, Solution<Scalar>*rsln, bool solutions_for_adapt,
      unsigned int error_flags)
    {
      this->tick();
      if(num != 1)
        throw Exceptions::LengthException(1, 1, num);
      double result = calc_err_internal(sln, rsln, NULL, solutions_for_adapt, error_flags);
      this->tick();
      this->info("Adaptivity: error estimate calculation duration: %f s.", this->last());
      return result;
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_est(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
      Hermes::vector<double>* component_errors, bool solutions_for_adapt,
      unsigned int error_flags)
    {
      this->tick();
      if(slns.size() != num)
        throw Exceptions::LengthException(1, slns.size(), num);
      if(rslns.size() != num)
        throw Exceptions::LengthException(2, rslns.size(), num);
      double result = calc_err_internal(slns, rslns, component_errors, solutions_for_adapt, error_flags);
      this->tick();
      this->info("Adaptivity: error estimate calculation duration: %f s.", this->last());
      return result;
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_exact(Solution<Scalar>*sln, Solution<Scalar>*rsln, bool solutions_for_adapt,
      unsigned int error_flags)
    {
      this->tick();
      if(num != 1)
        throw Exceptions::LengthException(1, 1, num);
      OGProjection<Scalar> ogProjection;
      typename Mesh::ReferenceMeshCreator ref_mesh_creator(this->spaces[0]->get_mesh());
      Mesh* ref_mesh = ref_mesh_creator.create_ref_mesh();
      typename Space<Scalar>::ReferenceSpaceCreator ref_space_creator(this->spaces[0], ref_mesh, 0);
      Space<Scalar>* ref_space = ref_space_creator.create_ref_space();
      Solution<Scalar> ref_sln_local;
      ogProjection.project_global(ref_space, rsln, &ref_sln_local);
      double result = calc_err_internal(sln, &ref_sln_local, NULL, solutions_for_adapt, error_flags);
      delete ref_space;
      delete ref_mesh;
      this->tick();
      this->info("Adaptivity: exact error calculation duration: %f s.", this->last());
      return result;
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_exact(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
      Hermes::vector<double>* component_errors, bool solutions_for_adapt,
      unsigned int error_flags)
    {
      this->tick();
      if(slns.size() != num)
        throw Exceptions::LengthException(1, slns.size(), num);
      if(rslns.size() != num)
        throw Exceptions::LengthException(2, rslns.size(), num);
      Mesh** ref_meshes = new Mesh*[num];
      Space<Scalar>** ref_spaces = new Space<Scalar>*[num];
      Hermes::vector<Solution<Scalar>*> ref_slns_local;
      for(unsigned int i = 0; i < num; i++)
      {
        OGProjection<Scalar> ogProjection;
        typename Mesh::ReferenceMeshCreator ref_mesh_creator(this->spaces[i]->get_mesh());
        ref_meshes[i] = ref_mesh_creator.create_ref_mesh();
        typename Space<Scalar>::ReferenceSpaceCreator ref_space_creator(this->spaces[i], ref_mesh_creator.create_ref_mesh(), 0);
        ref_spaces[i] = ref_space_creator.create_ref_space();
        ref_slns_local.push_back(new Solution<Scalar>);
        ogProjection.project_global(ref_space_creator.create_ref_space(), rslns[i], ref_slns_local.back());
      }
      double result = calc_err_internal(slns, ref_slns_local, component_errors, solutions_for_adapt, error_flags);
      for(unsigned int i = 0; i < num; i++)
      {
        delete ref_spaces[i];
        delete ref_meshes[i];
      }
      delete [] ref_meshes;
      delete [] ref_spaces;
      this->tick();
      this->info("Adaptivity: exact error calculation duration: %f s.", this->last());
      return result;
    }

    template<typename Scalar>
    bool Adapt<Scalar>::adapt(RefinementSelectors::Selector<Scalar>* refinement_selector, double thr, int strat,
      int regularize, double to_be_processed)
    {
      if(refinement_selector==NULL)
        throw Exceptions::NullException(1);
      Hermes::vector<RefinementSelectors::Selector<Scalar> *> refinement_selectors;
      refinement_selectors.push_back(refinement_selector);
      return adapt(refinement_selectors, thr, strat, regularize, to_be_processed);
    }

    template<typename Scalar>
    void Adapt<Scalar>::fix_shared_mesh_refinements(Mesh** meshes, std::vector<ElementToRefine>& elems_to_refine,
      int** idx, RefinementSelectors::Selector<Scalar> *** refinement_selectors)
    {
      int num_elem_to_proc = elems_to_refine.size();

      RefinementSelectors::Selector<Scalar>** current_refinement_selectors;

      for(int inx = 0; inx < num_elem_to_proc; inx++)
      {
        current_refinement_selectors = refinement_selectors[omp_get_thread_num()];
        ElementToRefine& elem_ref = elems_to_refine[inx];
        int current_quad_order = this->spaces[elem_ref.comp]->get_element_order(elem_ref.id);
        Element* current_elem = meshes[elem_ref.comp]->get_element(elem_ref.id);

        //select a refinement used by all components that share a mesh which is about to be refined
        int selected_refinement = elem_ref.split;
        for (int j = 0; j < this->num; j++)
        {
          if(selected_refinement == H2D_REFINEMENT_H) break; // iso refinement is max what can be recieved
          if(j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if a mesh is shared
            int ii = idx[elem_ref.id][j];
            if(ii >= 0) { // and the sample element is about to be refined by another compoment
              const ElementToRefine& elem_ref_ii = elems_to_refine[ii];
              if(elem_ref_ii.split != selected_refinement && elem_ref_ii.split != H2D_REFINEMENT_P) { //select more complicated refinement
                if((elem_ref_ii.split == H2D_REFINEMENT_ANISO_H || elem_ref_ii.split == H2D_REFINEMENT_ANISO_V) && selected_refinement == H2D_REFINEMENT_P)
                  selected_refinement = elem_ref_ii.split;
                else
                  selected_refinement = H2D_REFINEMENT_H;
              }
            }
          }
        }

        //fix other refinements according to the selected refinement
        if(selected_refinement != H2D_REFINEMENT_P)
        {
          //get suggested orders for the selected refinement
          const int* suggested_orders = NULL;
          if(selected_refinement == H2D_REFINEMENT_H)
            suggested_orders = elem_ref.q;

          //update orders
          for (int j = 0; j < this->num; j++)
          {
            if(j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if components share the mesh
              // change currently processed refinement
              if(elem_ref.split != selected_refinement)
              {
                elem_ref.split = selected_refinement;
                current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref.split, elem_ref.p, suggested_orders);
              }

              // change other refinements
              int ii = idx[elem_ref.id][j];
              if(ii >= 0)
              {
                ElementToRefine& elem_ref_ii = elems_to_refine[ii];
                if(elem_ref_ii.split != selected_refinement)
                {
                  elem_ref_ii.split = selected_refinement;
                  current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_ii.split, elem_ref_ii.p, suggested_orders);
                }
              }
              else
              { // element (of the other comp.) not refined at all: assign refinement
                ElementToRefine elem_ref_new(elem_ref.id, j);
                elem_ref_new.split = selected_refinement;
                current_refinement_selectors[j]->generate_shared_mesh_orders(current_elem, current_quad_order, elem_ref_new.split, elem_ref_new.p, suggested_orders);
                elems_to_refine.push_back(elem_ref_new);
              }
            }
          }
        }
      }
    }

    template<typename Scalar>
    double Adapt<Scalar>::get_element_error_squared(int component, int id) const
    {
      if(!have_errors)
        throw Exceptions::Exception("element errors have to be calculated first, call Adapt<Scalar>::calc_err_est().");
      return errors[component][id];
    };

    template<typename Scalar>
    const Hermes::vector<typename Adapt<Scalar>::ElementReference>& Adapt<Scalar>::get_regular_queue() const
    {
      return regular_queue;
    };

    template<typename Scalar>
    void Adapt<Scalar>::homogenize_shared_mesh_orders(Mesh** meshes)
    {
      Element* e;
      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, meshes[i])
        {
          int current_quad_order = this->spaces[i]->get_element_order(e->id);
          int current_order_h = H2D_GET_H_ORDER(current_quad_order), current_order_v = H2D_GET_V_ORDER(current_quad_order);

          for (int j = 0; j < this->num; j++)
            if((j != i) && (meshes[j] == meshes[i])) // components share the mesh
            {
              int quad_order = this->spaces[j]->get_element_order(e->id);
              current_order_h = std::max(current_order_h, H2D_GET_H_ORDER(quad_order));
              current_order_v = std::max(current_order_v, H2D_GET_V_ORDER(quad_order));
            }

            this->spaces[i]->set_element_order_internal(e->id, H2D_MAKE_QUAD_ORDER(current_order_h, current_order_v));
        }
      }
    }

    template<typename Scalar>
    const std::vector<ElementToRefine>& Adapt<Scalar>::get_last_refinements() const
    {
      return last_refinements;
    }

    template<typename Scalar>
    void Adapt<Scalar>::apply_refinements(std::vector<ElementToRefine>& elems_to_refine)
    {
      for (std::vector<ElementToRefine>::const_iterator elem_ref = elems_to_refine.begin();
        elem_ref != elems_to_refine.end(); elem_ref++)
        apply_refinement(*elem_ref);
    }

    template<typename Scalar>
    void Adapt<Scalar>::apply_refinement(const ElementToRefine& elem_ref)
    {
      Space<Scalar>* space = this->spaces[elem_ref.comp];
      Mesh* mesh = space->get_mesh();

      Element* e;
      e = mesh->get_element(elem_ref.id);

      if(elem_ref.split == H2D_REFINEMENT_P)
      {
        space->set_element_order_internal(elem_ref.id, elem_ref.p[0]);
        space->edata[elem_ref.id].changed_in_last_adaptation = true;
      }
      else if(elem_ref.split == H2D_REFINEMENT_H)
      {
        if(e->active)
          mesh->refine_element_id(elem_ref.id);
        for (int j = 0; j < 4; j++)
        {
          space->set_element_order_internal(e->sons[j]->id, elem_ref.p[j]);
          space->edata[e->sons[j]->id].changed_in_last_adaptation = true;
        }
      }
      else
      {
        if(e->active)
          mesh->refine_element_id(elem_ref.id, elem_ref.split);
        for (int j = 0; j < 2; j++)
        {
          space->set_element_order_internal(e->sons[ (elem_ref.split == 1) ? j : j + 2 ]->id, elem_ref.p[j]);
          space->edata[e->sons[ (elem_ref.split == 1) ? j : j + 2 ]->id].changed_in_last_adaptation = true;
        }
      }
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_error_form(int i, int j, typename Adapt<Scalar>::MatrixFormVolError* form)
    {
      if(form->i < 0 || form->i >= this->num)
        throw Exceptions::ValueException("component number", form->i, 0, this->num);
      if(form->j < 0 || form->j >= this->num)
        throw Exceptions::ValueException("component number", form->j, 0, this->num);

      // FIXME: Memory leak - always for i == j (see the constructor), may happen for i != j
      //        if user does not delete previously set error forms by himself.
      if(own_forms[i][j] && error_form[i][j] != NULL)
        delete error_form[i][j];
      error_form[i][j] = form;
      norm_form[i][j] = error_form[i][j];
      own_forms[i][j] = false;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_error_form(typename Adapt<Scalar>::MatrixFormVolError* form)
    {
      set_error_form(0, 0, form);
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_norm_form(int i, int j, typename Adapt<Scalar>::MatrixFormVolError* form)
    {
      if(form->i < 0 || form->i >= this->num)
        throw Exceptions::ValueException("component number", form->i, 0, this->num);
      if(form->j < 0 || form->j >= this->num)
        throw Exceptions::ValueException("component number", form->j, 0, this->num);

      norm_form[i][j] = form;
    }

    template<typename Scalar>
    void Adapt<Scalar>::set_norm_form(typename Adapt<Scalar>::MatrixFormVolError* form)
    {
      set_norm_form(0, 0, form);
    }

    template<typename Scalar>
    double Adapt<Scalar>::eval_error(typename Adapt<Scalar>::MatrixFormVolError* form,
      MeshFunction<Scalar>*sln1, MeshFunction<Scalar>*sln2, MeshFunction<Scalar>*rsln1,
      MeshFunction<Scalar>*rsln2)
    {
      RefMap *rv1 = sln1->get_refmap();
      RefMap *rv2 = sln2->get_refmap();
      RefMap *rrv1 = rsln1->get_refmap();
      RefMap *rrv2 = rsln2->get_refmap();

      // determine the integration order
      int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
      Func<Hermes::Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
      Func<Hermes::Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

      double fake_wt = 1.0;
      Geom<Hermes::Ord>* fake_e = init_geom_ord();
      Hermes::Ord o = form->ord(1, &fake_wt, NULL, ou, ov, fake_e, NULL);
      int order = rrv1->get_inv_ref_order();
      order += o.get_order();
      if(static_cast<Solution<Scalar>*>(rsln1) || static_cast<Solution<Scalar>*>(rsln2))
      {
        if(static_cast<Solution<Scalar>*>(rsln1)->get_type() == HERMES_EXACT)
        { limit_order_nowarn(order, rv1->get_active_element()->get_mode()); }
        else
          limit_order(order, rv1->get_active_element()->get_mode());
      }
      else
        limit_order(order, rv1->get_active_element()->get_mode());

      ou->free_ord(); delete ou;
      ov->free_ord(); delete ov;
      delete fake_e;

      // eval the form
      Quad2D* quad = sln1->get_quad_2d();
      double3* pt = quad->get_points(order, sln1->get_active_element()->get_mode());
      int np = quad->get_num_points(order, sln1->get_active_element()->get_mode());

      // init geometry and jacobian*weights
      Geom<double>* e = init_geom_vol(rrv1, order);
      double* jac = rrv1->get_jacobian(order);
      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * jac[i];

      // function values and values of external functions
      Func<Scalar>* err1 = init_fn(sln1, order);
      Func<Scalar>* err2 = init_fn(sln2, order);
      Func<Scalar>* v1 = init_fn(rsln1, order);
      Func<Scalar>* v2 = init_fn(rsln2, order);

      err1->subtract(v1);
      err2->subtract(v2);

      Scalar res = form->value(np, jwt, NULL, err1, err2, e, NULL);

      e->free(); delete e;
      delete [] jwt;
      err1->free_fn(); delete err1;
      err2->free_fn(); delete err2;
      v1->free_fn(); delete v1;
      v2->free_fn(); delete v2;

      return std::abs(res);
    }

    template<typename Scalar>
    double Adapt<Scalar>::eval_error_norm(typename Adapt<Scalar>::MatrixFormVolError* form,
      MeshFunction<Scalar>*rsln1, MeshFunction<Scalar>*rsln2)
    {
      RefMap *rrv1 = rsln1->get_refmap();
      RefMap *rrv2 = rsln2->get_refmap();

      // determine the integration order
      int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
      Func<Hermes::Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
      Func<Hermes::Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

      double fake_wt = 1.0;
      Geom<Hermes::Ord>* fake_e = init_geom_ord();
      Hermes::Ord o = form->ord(1, &fake_wt, NULL, ou, ov, fake_e, NULL);
      int order = rrv1->get_inv_ref_order();
      order += o.get_order();
      if(static_cast<Solution<Scalar>*>(rsln1) || static_cast<Solution<Scalar>*>(rsln2))
      {
        if(static_cast<Solution<Scalar>*>(rsln1)->get_type() == HERMES_EXACT)
        { limit_order_nowarn(order, rrv1->get_active_element()->get_mode());  }
        else
          limit_order(order, rrv1->get_active_element()->get_mode());
      }
      else
        limit_order(order, rrv1->get_active_element()->get_mode());

      ou->free_ord(); delete ou;
      ov->free_ord(); delete ov;
      delete fake_e;

      // eval the form
      Quad2D* quad = rsln1->get_quad_2d();
      double3* pt = quad->get_points(order, rrv1->get_active_element()->get_mode());
      int np = quad->get_num_points(order, rrv1->get_active_element()->get_mode());

      // init geometry and jacobian*weights
      Geom<double>* e = init_geom_vol(rrv1, order);
      double* jac = rrv1->get_jacobian(order);
      double* jwt = new double[np];
      for(int i = 0; i < np; i++)
        jwt[i] = pt[i][2] * jac[i];

      // function values
      Func<Scalar>* v1 = init_fn(rsln1, order);
      Func<Scalar>* v2 = init_fn(rsln2, order);

      Scalar res = form->value(np, jwt, NULL, v1, v2, e, NULL);

      e->free(); delete e;
      delete [] jwt;
      v1->free_fn(); delete v1;
      v2->free_fn(); delete v2;

      return std::abs(res);
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_internal(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
      Hermes::vector<double>* component_errors, bool solutions_for_adapt, unsigned int error_flags)
    {
      int i, j;
      
      bool compatible_meshes = true;
      for (int space_i = 0; space_i < this->num; space_i++)
      {
        Element* e;
        for_all_active_elements(e, slns[space_i]->get_mesh())
        {
          Element* e_ref = rslns[space_i]->get_mesh()->get_element(e->id);
          if(e_ref == NULL)
          {
            compatible_meshes = false;
            break;
          }
          if(!e_ref->active)
          {
            if(e_ref->sons[0] == NULL || e_ref->sons[2] == NULL)
            {
              compatible_meshes = false;
              break;
            }
            if(!e_ref->sons[0]->active || !e_ref->sons[2]->active)
            {
              compatible_meshes = false;
              break;
            }
          }
        }
      }

      if(!compatible_meshes)
        throw Exceptions::Exception("Reference space not created by an isotropic (p-, h-, or hp-) refinement from the coarse space.");

      if(slns.size() != this->num)
        throw Exceptions::LengthException(0, slns.size(), this->num);

      Solution<Scalar>* rslns_original[H2D_MAX_COMPONENTS];
      Solution<Scalar>* slns_original[H2D_MAX_COMPONENTS];

      for (i = 0; i < this->num; i++)
      {
        slns_original[i] = this->sln[i];
        this->sln[i] = slns[i];
        sln[i]->set_quad_2d(&g_quad_2d_std);
      }
      for (i = 0; i < this->num; i++)
      {
        rslns_original[i] = this->rsln[i];
        this->rsln[i] = rslns[i];
        rsln[i]->set_quad_2d(&g_quad_2d_std);
      }

      have_coarse_solutions = true;
      have_reference_solutions = true;

      // Prepare multi-mesh traversal and error arrays.
      const Mesh **meshes = new const Mesh *[2 * num];
      Transformable **tr = new Transformable *[2 * num];
      Traverse trav(true);
      num_act_elems = 0;
      for (i = 0; i < num; i++)
      {
        meshes[i] = sln[i]->get_mesh();
        meshes[i + num] = rsln[i]->get_mesh();
        tr[i] = sln[i];
        tr[i + num] = rsln[i];

        num_act_elems += sln[i]->get_mesh()->get_num_active_elements();

        int max = meshes[i]->get_max_element_id();
        if(solutions_for_adapt)
        {
          if(errors[i] != NULL) delete [] errors[i];
          errors[i] = new double[max];
          memset(errors[i], 0, sizeof(double) * max);
        }
      }

      double total_norm = 0.0;
      double *norms = new double[num];
      memset(norms, 0, num * sizeof(double));
      double *errors_components = new double[num];
      memset(errors_components, 0, num * sizeof(double));
      if(solutions_for_adapt) this->errors_squared_sum = 0.0;
      double total_error = 0.0;

      // Calculate error.
      Traverse::State * ee;
      trav.begin(2 * num, meshes, tr);
      while ((ee = trav.get_next_state()) != NULL)
      {
        for (i = 0; i < num; i++)
        {
          for (j = 0; j < num; j++)
          {
            if(error_form[i][j] != NULL)
            {
              double err, nrm;
              err = eval_error(error_form[i][j], sln[i], sln[j], rsln[i], rsln[j]);
              nrm = eval_error_norm(norm_form[i][j], rsln[i], rsln[j]);

              norms[i] += nrm;
              total_norm  += nrm;
              total_error += err;
              errors_components[i] += err;
              if(solutions_for_adapt)
                this->errors[i][ee->e[i]->id] += err;
            }
          }
        }
      }
      trav.finish();

      // Store the calculation for each solution component separately.
      if(component_errors != NULL)
      {
        component_errors->clear();
        for (int i = 0; i < num; i++)
        {
          if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
            component_errors->push_back(sqrt(errors_components[i]));
          else if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
            component_errors->push_back(sqrt(errors_components[i]/norms[i]));
          else
          {
            throw Hermes::Exceptions::Exception("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
            return -1.0;
          }
        }
      }

      // Make the error relative if needed.
      if(solutions_for_adapt)
      {
        if((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
        {
          for (int i = 0; i < this->num; i++)
          {
            Element* e;
            for_all_active_elements(e, meshes[i])
            {
              errors[i][e->id] /= norms[i];
            }
          }
        }

        this->errors_squared_sum = total_error;

        // Element error mask is used here, because this variable is used in the adapt()
        // function, where the processed error (sum of errors of processed element errors)
        // is matched to this variable.
        if((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
          errors_squared_sum = errors_squared_sum / total_norm;
      }

      // Prepare an ordered list of elements according to an error.
      if(solutions_for_adapt)
      {
        fill_regular_queue(meshes);
        have_errors = true;
      }
      else
      {
        for (i = 0; i < this->num; i++)
        {
          this->sln[i] = slns_original[i];
          this->rsln[i] = rslns_original[i];
        }
      }

      delete [] meshes;
      delete [] tr;
      delete [] norms;
      delete [] errors_components;

      // Return error value.
      if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
        return sqrt(total_error);
      else if((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
        return sqrt(total_error / total_norm);
      else
      {
        throw Hermes::Exceptions::Exception("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
        return -1.0;
      }
    }

    template<typename Scalar>
    double Adapt<Scalar>::calc_err_internal(Solution<Scalar>* sln, Solution<Scalar>* rsln,
      Hermes::vector<double>* component_errors, bool solutions_for_adapt,
      unsigned int error_flags)
    {
      Hermes::vector<Solution<Scalar>*> slns;
      slns.push_back(sln);
      Hermes::vector<Solution<Scalar>*> rslns;
      rslns.push_back(rsln);
      return calc_err_internal(slns, rslns, component_errors, solutions_for_adapt, error_flags);
    }

    template<typename Scalar>
    Adapt<Scalar>::CompareElements::CompareElements(double** errors): errors(errors)
    {
    }

    template<typename Scalar>
    bool Adapt<Scalar>::CompareElements::operator()(const ElementReference& e1, const ElementReference& e2) const
    {
      return errors[e1.comp][e1.id] > errors[e2.comp][e2.id];
    }

    template<typename Scalar>
    void Adapt<Scalar>::fill_regular_queue(const Mesh** meshes)
    {
      //prepare space for queue (it is assumed that it will only grow since we can just split)
      regular_queue.clear();
      if(num_act_elems < (int)regular_queue.capacity())
      {
        Hermes::vector<ElementReference> empty_refs;
        regular_queue.swap(empty_refs); //deallocate
        regular_queue.reserve(num_act_elems); //allocate
      }

      //prepare initial fill
      Element* e;
      typename Hermes::vector<ElementReference>::iterator elem_info = regular_queue.begin();
      for (int i = 0; i < this->num; i++)
      {
        for_all_active_elements(e, meshes[i])
        {
          regular_queue.push_back(ElementReference(e->id, i));
        }
      }
      //sort
      std::sort(regular_queue.begin(), regular_queue.end(), CompareElements(errors));
    }

    template HERMES_API class Adapt<double>;
    template HERMES_API class Adapt<std::complex<double> >;
  }
}
