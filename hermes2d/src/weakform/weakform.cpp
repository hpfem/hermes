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

#include "global.h"
#include "api2d.h"
#include "weakform.h"
#include "matrix.h"
#include "forms.h"
#include "shapeset_l2_all.h"
#include "shapeset_hc_all.h"
#include "shapeset_hd_all.h"
#include "shapeset_h1_all.h"
using namespace Hermes::Algebra::DenseMatrixOperations;
namespace Hermes
{
  namespace Hermes2D
  {

    /// This is to be used by weak forms specifying numerical flux through interior edges.
    /// Forms with this identifier will receive DiscontinuousFunc representations of shape
    /// and ext. functions, which they may query for values on either side of given interface.
    static const std::string H2D_DG_INNER_EDGE = "-1234567";

    template<typename Scalar>
    WeakForm<Scalar>::WeakForm(unsigned int neq, bool mat_free) : Hermes::Mixins::Loggable(true), warned_nonOverride(false)
    {
      this->neq = neq;
      this->is_matfree = mat_free;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::free_ext()
    {
      for(unsigned int i = 0; i < this->ext.size(); i++)
        delete this->ext[i];
      for(unsigned int i = 0; i < this->forms.size(); i++)
        for(unsigned int j = 0; j < get_forms()[i]->ext.size(); j++)
          delete get_forms()[i]->ext[j];
    }

    template<typename Scalar>
    WeakForm<Scalar>::~WeakForm()
    {
      for(unsigned int i = 0; i < this->forms.size(); i++)
        delete get_forms()[i];
      delete_all();
    }

    template<typename Scalar>
    WeakForm<Scalar>* WeakForm<Scalar>::clone() const
    {
      if(!this->warned_nonOverride)
#pragma omp critical (warning_weakform_nonOverride)
      {
        if(!this->warned_nonOverride)
          this->warn("Using default WeakForm<Scalar>::clone, if you have any dynamically created data in your WeakForm constructor, you need to overload this method!");
        const_cast<WeakForm<Scalar>*>(this)->warned_nonOverride = true;
      }
      return new WeakForm(*this);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_active_state(Element** e)
    {
      // Nothing here, supposed to be overriden if necessary.
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_active_edge_state(Element** e, int isurf)
    {
      // Nothing here, supposed to be overriden if necessary.
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_active_DG_state(Element** e, int isurf)
    {
      // Nothing here, supposed to be overriden if necessary.
    }

    template<typename Scalar>
    void WeakForm<Scalar>::cloneMembers(const WeakForm<Scalar>* otherWf)
    {
      this->mfvol.clear();
      this->mfsurf.clear();
      this->mfDG.clear();
      this->vfvol.clear();
      this->vfsurf.clear();
      this->vfDG.clear();
      this->forms.clear();
      this->ext.clear();

      for(unsigned int i = 0; i < otherWf->forms.size(); i++)
      {
        if(dynamic_cast<MatrixFormVol<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<MatrixFormVol<Scalar>*>(otherWf->forms[i]))->clone());
        if(dynamic_cast<MatrixFormSurf<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<MatrixFormSurf<Scalar>*>(otherWf->forms[i]))->clone());
        if(dynamic_cast<MatrixFormDG<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<MatrixFormDG<Scalar>*>(otherWf->forms[i]))->clone());

        if(dynamic_cast<VectorFormVol<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<VectorFormVol<Scalar>*>(otherWf->forms[i]))->clone());
        if(dynamic_cast<VectorFormSurf<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<VectorFormSurf<Scalar>*>(otherWf->forms[i]))->clone());
        if(dynamic_cast<VectorFormDG<Scalar>*>(otherWf->forms[i]) != NULL)
          this->forms.push_back((dynamic_cast<VectorFormDG<Scalar>*>(otherWf->forms[i]))->clone());

        Hermes::vector<MeshFunction<Scalar>*> newExt;
        for(unsigned int ext_i = 0; ext_i < otherWf->forms[i]->ext.size(); ext_i++)
          newExt.push_back(otherWf->forms[i]->ext[ext_i]->clone());
        this->forms.back()->set_ext(newExt);
        this->forms.back()->wf = this;

        if(dynamic_cast<MatrixFormVol<Scalar>*>(otherWf->forms[i]) != NULL)
          this->mfvol.push_back(dynamic_cast<MatrixFormVol<Scalar>*>(this->forms.back()));
        if(dynamic_cast<MatrixFormSurf<Scalar>*>(otherWf->forms[i]) != NULL)
          this->mfsurf.push_back(dynamic_cast<MatrixFormSurf<Scalar>*>(this->forms.back()));
        if(dynamic_cast<MatrixFormDG<Scalar>*>(otherWf->forms[i]) != NULL)
          this->mfDG.push_back(dynamic_cast<MatrixFormDG<Scalar>*>(this->forms.back()));

        if(dynamic_cast<VectorFormVol<Scalar>*>(otherWf->forms[i]) != NULL)
          this->vfvol.push_back(dynamic_cast<VectorFormVol<Scalar>*>(this->forms.back()));
        if(dynamic_cast<VectorFormSurf<Scalar>*>(otherWf->forms[i]) != NULL)
          this->vfsurf.push_back(dynamic_cast<VectorFormSurf<Scalar>*>(this->forms.back()));
        if(dynamic_cast<VectorFormDG<Scalar>*>(otherWf->forms[i]) != NULL)
          this->vfDG.push_back(dynamic_cast<VectorFormDG<Scalar>*>(this->forms.back()));
      }
      for(unsigned int i = 0; i < otherWf->ext.size(); i++)
        this->ext.push_back(otherWf->ext[i]->clone());
    }

    template<typename Scalar>
    void WeakForm<Scalar>::delete_all()
    {
      mfvol.clear();
      mfsurf.clear();
      mfDG.clear();
      vfvol.clear();
      vfsurf.clear();
      vfDG.clear();
      forms.clear();
    };

    template<typename Scalar>
    void WeakForm<Scalar>::set_ext(MeshFunction<Scalar>* ext)
    {
      this->ext.clear();
      this->ext.push_back(ext);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_ext(Hermes::vector<MeshFunction<Scalar>*> ext)
    {
      this->ext = ext;
    }

    template<typename Scalar>
    Hermes::vector<MeshFunction<Scalar>*> WeakForm<Scalar>::get_ext() const
    {
      return this->ext;
    }

    template<typename Scalar>
    Form<Scalar>::Form() : scaling_factor(1.0), u_ext_offset(0), wf(NULL)
    {
      areas.push_back(HERMES_ANY);
      stage_time = 0.0;
    }

    template<typename Scalar>
    Form<Scalar>::~Form()
    {
    }

    template<typename Scalar>
    void Form<Scalar>::set_current_stage_time(double time)
    {
      stage_time = time;
    }

    template<typename Scalar>
    double Form<Scalar>::get_current_stage_time() const
    {
      return stage_time;
    }

    template<typename Scalar>
    void Form<Scalar>::set_area(std::string area)
    {
      areas.clear();
      areas.push_back(area);
    }
    template<typename Scalar>
    void Form<Scalar>::set_areas(Hermes::vector<std::string> areas)
    {
      this->areas = areas;
    }

    template<typename Scalar>
    Hermes::vector<std::string> Form<Scalar>::getAreas() const
    {
      return this->areas;
    }

    template<typename Scalar>
    void Form<Scalar>::setScalingFactor(double scalingFactor)
    {
      this->scaling_factor = scalingFactor;
    }

    template<typename Scalar>
    void Form<Scalar>::set_uExtOffset(int u_ext_offset)
    {
      this->u_ext_offset = u_ext_offset;
    }

    template<typename Scalar>
    void Form<Scalar>::set_ext(MeshFunction<Scalar>* ext)
    {
      this->ext.clear();
      this->ext.push_back(ext);
    }

    template<typename Scalar>
    void Form<Scalar>::set_ext(Hermes::vector<MeshFunction<Scalar>*> ext)
    {
      this->ext = ext;
    }

    template<typename Scalar>
    Hermes::vector<MeshFunction<Scalar>*> Form<Scalar>::get_ext() const
    {
      return this->ext;
    }

    template<typename Scalar>
    MatrixForm<Scalar>::MatrixForm(unsigned int i, unsigned int j) :
    Form<Scalar>(), sym(HERMES_NONSYM), i(i), j(j), previous_iteration_space_index(-1)
    {
    }

    template<typename Scalar>
    MatrixForm<Scalar>::~MatrixForm()
    {
    }

    template<typename Scalar>
    Scalar MatrixForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, Func<Scalar> **ext) const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("MatrixForm<Scalar>::value");
      return 0.0;
    }

    template<typename Scalar>
    Hermes::Ord MatrixForm<Scalar>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, Func<Ord> **ext) const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("MatrixForm<Scalar>::ord");
      return Hermes::Ord();
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>::MatrixFormVol(unsigned int i, unsigned int j) :
    MatrixForm<Scalar>(i, j)
    {
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>::~MatrixFormVol()
    {
    }

    template<typename Scalar>
    void MatrixFormVol<Scalar>::setSymFlag(SymFlag sym)
    {
      this->sym = sym;
    }

    template<typename Scalar>
    SymFlag MatrixFormVol<Scalar>::getSymFlag() const
    {
      return this->sym;
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>* MatrixFormVol<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("MatrixFormVol<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>::MatrixFormSurf(unsigned int i, unsigned int j) :
    MatrixForm<Scalar>(i, j)
    {
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>* MatrixFormSurf<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("MatrixFormSurf<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>::~MatrixFormSurf()
    {
    }

    template<typename Scalar>
    MatrixFormDG<Scalar>::MatrixFormDG(unsigned int i, unsigned int j) :
    MatrixForm<Scalar>(i, j)
    {
      this->set_area(H2D_DG_INNER_EDGE);
    }

    template<typename Scalar>
    MatrixFormDG<Scalar>::~MatrixFormDG()
    {
    }

    template<typename Scalar>
    MatrixFormDG<Scalar>* MatrixFormDG<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("MatrixFormDG<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    VectorForm<Scalar>::VectorForm(unsigned int i) :
    Form<Scalar>(), i(i)
    {
    }

    template<typename Scalar>
    VectorForm<Scalar>::~VectorForm()
    {
    }

    template<typename Scalar>
    VectorFormVol<Scalar>::VectorFormVol(unsigned int i) :
    VectorForm<Scalar>(i)
    {
    }

    template<typename Scalar>
    VectorFormVol<Scalar>::~VectorFormVol()
    {
    }

    template<typename Scalar>
    Scalar VectorForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
      Geom<double> *e, Func<Scalar> **ext) const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("VectorForm<Scalar>::value");
      return 0.0;
    }

    template<typename Scalar>
    Hermes::Ord VectorForm<Scalar>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, Func<Ord> **ext) const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("VectorForm<Scalar>::ord");
      return Hermes::Ord();
    }

    template<typename Scalar>
    VectorFormVol<Scalar>* VectorFormVol<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("VectorFormVol<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    VectorFormSurf<Scalar>::VectorFormSurf(unsigned int i) :
    VectorForm<Scalar>(i)
    {
    }

    template<typename Scalar>
    VectorFormSurf<Scalar>::~VectorFormSurf()
    {
    }

    template<typename Scalar>
    VectorFormSurf<Scalar>* VectorFormSurf<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("VectorFormSurf<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    VectorFormDG<Scalar>::VectorFormDG(unsigned int i) :
    VectorForm<Scalar>(i)
    {
      this->set_area(H2D_DG_INNER_EDGE);
    }

    template<typename Scalar>
    VectorFormDG<Scalar>::~VectorFormDG()
    {
    }

    template<typename Scalar>
    VectorFormDG<Scalar>* VectorFormDG<Scalar>::clone() const
    {
      throw Hermes::Exceptions::MethodNotOverridenException("VectorFormDG<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_matrix_form(MatrixFormVol<Scalar>* form)
    {
      if(form->i >= neq || form->j >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");
      if(form->sym < -1 || form->sym > 1)
        throw Hermes::Exceptions::Exception("\"sym\" must be -1, 0 or 1.");
      if(form->sym < 0 && form->i == form->j)
        throw Hermes::Exceptions::Exception("Only off-diagonal forms can be antisymmetric.");
      if(mfvol.size() > 100)
      {
        this->warn("Large number of forms (> 100). Is this the intent?");
      }

      form->set_weakform(this);
      mfvol.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_matrix_form_surf(MatrixFormSurf<Scalar>* form)
    {
      if(form->i >= neq || form->j >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      mfsurf.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_matrix_form_DG(MatrixFormDG<Scalar>* form)
    {
      if(form->i >= neq || form->j >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      mfDG.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_vector_form(VectorFormVol<Scalar>* form)
    {
      if(form->i >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");
      form->set_weakform(this);
      vfvol.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_vector_form_surf(VectorFormSurf<Scalar>* form)
    {
      if(form->i >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      vfsurf.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_vector_form_DG(VectorFormDG<Scalar>* form)
    {
      if(form->i >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      vfDG.push_back(form);
      forms.push_back(form);
    }

    template<typename Scalar>
    Hermes::vector<Form<Scalar> *> WeakForm<Scalar>::get_forms() const
    {
      return forms;
    }

    template<typename Scalar>
    Hermes::vector<MatrixFormVol<Scalar> *> WeakForm<Scalar>::get_mfvol() const
    {
      return mfvol;
    }
    template<typename Scalar>
    Hermes::vector<MatrixFormSurf<Scalar> *> WeakForm<Scalar>::get_mfsurf() const
    {
      return mfsurf;
    }
    template<typename Scalar>
    Hermes::vector<MatrixFormDG<Scalar> *> WeakForm<Scalar>::get_mfDG() const
    {
      return mfDG;
    }
    template<typename Scalar>
    Hermes::vector<VectorFormVol<Scalar> *> WeakForm<Scalar>::get_vfvol() const
    {
      return vfvol;
    }
    template<typename Scalar>
    Hermes::vector<VectorFormSurf<Scalar> *> WeakForm<Scalar>::get_vfsurf() const
    {
      return vfsurf;
    }
    template<typename Scalar>
    Hermes::vector<VectorFormDG<Scalar> *> WeakForm<Scalar>::get_vfDG() const
    {
      return vfDG;
    }

    template<typename Scalar>
    bool** WeakForm<Scalar>::get_blocks(bool force_diagonal_blocks) const
    {
      bool** blocks = new_matrix<bool>(neq, neq);
      for (unsigned int i = 0; i < neq; i++)
      {
        for (unsigned int j = 0; j < neq; j++)
          blocks[i][j] = false;
        if(force_diagonal_blocks)
          blocks[i][i] = true;
      }
      for (unsigned i = 0; i < mfvol.size(); i++)
      {
        if(fabs(mfvol[i]->scaling_factor) > 1e-12)
          blocks[mfvol[i]->i][mfvol[i]->j] = true;
        if(mfvol[i]->sym)
          if(fabs(mfvol[i]->scaling_factor) > 1e-12)
            blocks[mfvol[i]->j][mfvol[i]->i] = true;
      }
      for (unsigned i = 0; i < mfsurf.size(); i++)
      {
        if(fabs(mfsurf[i]->scaling_factor) > 1e-12)
          blocks[mfsurf[i]->i][mfsurf[i]->j] = true;
      }

      return blocks;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_current_time(double time)
    {
      current_time = time;
    }

    template<typename Scalar>
    double WeakForm<Scalar>::get_current_time() const
    {
      return current_time;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::set_current_time_step(double time_step)
    {
      current_time_step = time_step;
    }

    template<typename Scalar>
    double WeakForm<Scalar>::get_current_time_step() const
    {
      return current_time_step;
    }

    template class HERMES_API WeakForm<double>;
    template class HERMES_API WeakForm<std::complex<double> >;
    template class HERMES_API Form<double>;
    template class HERMES_API Form<std::complex<double> >;
    template class HERMES_API MatrixForm<double>;
    template class HERMES_API MatrixForm<std::complex<double> >;
    template class HERMES_API MatrixFormVol<double>;
    template class HERMES_API MatrixFormVol<std::complex<double> >;
    template class HERMES_API MatrixFormSurf<double>;
    template class HERMES_API MatrixFormSurf<std::complex<double> >;
    template class HERMES_API MatrixFormDG<double>;
    template class HERMES_API MatrixFormDG<std::complex<double> >;
    template class HERMES_API VectorForm<double>;
    template class HERMES_API VectorForm<std::complex<double> >;
    template class HERMES_API VectorFormVol<double>;
    template class HERMES_API VectorFormVol<std::complex<double> >;
    template class HERMES_API VectorFormSurf<double>;
    template class HERMES_API VectorFormSurf<std::complex<double> >;
    template class HERMES_API VectorFormDG<double>;
    template class HERMES_API VectorFormDG<std::complex<double> >;
  }
}
