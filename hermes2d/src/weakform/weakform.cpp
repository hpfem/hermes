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
#include "weakform.h"
#include "matrix.h"
#include "forms.h"
#include "space.h"

using namespace Hermes::Algebra::DenseMatrixOperations;
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    WeakForm<Scalar>::WeakForm(unsigned int neq, bool mat_free)
    {
      this->neq = neq;
      this->seq = 0;
      this->is_matfree = mat_free;
    }

    template<typename Scalar>
    WeakForm<Scalar>::~WeakForm()
    {
      delete_all();
    }

    template<typename Scalar>
    void WeakForm<Scalar>::delete_all()
    {
      mfvol.clear();
      mfsurf.clear();
      vfvol.clear();
      vfsurf.clear();
    };

    template<typename Scalar>
    Form<Scalar>::Form(std::string area, Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    ext(ext), scaling_factor(scaling_factor), u_ext_offset(u_ext_offset)
    {
      areas.push_back(area);
      stage_time = 0.0;
    }

    template<typename Scalar>
    Form<Scalar>::Form(Hermes::vector<std::string> areas, Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    ext(ext), scaling_factor(scaling_factor), u_ext_offset(u_ext_offset)
    {
      this->areas = areas;
      stage_time = 0.0;
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
    MatrixForm<Scalar>::MatrixForm(unsigned int i, unsigned int j,
      std::string area, Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    Form<Scalar>(area, ext, scaling_factor, u_ext_offset), sym(0)
    {
      this->i = i;
      this->j = j;
    }

    template<typename Scalar>
    MatrixForm<Scalar>::MatrixForm(unsigned int i, unsigned int j,
      Hermes::vector<std::string> areas, Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    Form<Scalar>(areas, ext, scaling_factor, u_ext_offset), sym(0)
    {
      this->i = i;
      this->j = j;
    }

    template<typename Scalar>
    Scalar MatrixForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<Scalar> *ext) const
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("MatrixForm<Scalar>::value");
      return 0.0;
    }

    template<typename Scalar>
    Hermes::Ord MatrixForm<Scalar>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("MatrixForm<Scalar>::ord");
      return Hermes::Ord();
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>::MatrixFormVol(unsigned int i, unsigned int j,
      std::string area, SymFlag sym, Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    MatrixForm<Scalar>(i, j, area, ext, scaling_factor, u_ext_offset)
    {
      this->sym = sym;
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>::MatrixFormVol(unsigned int i, unsigned int j,
      Hermes::vector<std::string> areas, SymFlag sym, Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    MatrixForm<Scalar>(i, j, areas, ext, scaling_factor, u_ext_offset)
    {
      this->sym = sym;
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>* MatrixFormVol<Scalar>::clone()
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("MatrixFormVol<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>::MatrixFormSurf(unsigned int i, unsigned int j, std::string area,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    MatrixForm<Scalar>(i, j, area, ext, scaling_factor, u_ext_offset)
    {
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>::MatrixFormSurf(unsigned int i, unsigned int j, Hermes::vector<std::string> areas,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    MatrixForm<Scalar>(i, j, areas, ext, scaling_factor, u_ext_offset)
    {
    }

    template<typename Scalar>
    MatrixFormSurf<Scalar>* MatrixFormSurf<Scalar>::clone()
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("MatrixFormSurf<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    VectorForm<Scalar>::VectorForm(unsigned int i, std::string area,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    Form<Scalar>(area, ext, scaling_factor, u_ext_offset)
    {
      this->i = i;
    }

    template<typename Scalar>
    VectorForm<Scalar>::VectorForm(unsigned int i, Hermes::vector<std::string> areas,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    Form<Scalar>(areas, ext, scaling_factor, u_ext_offset)
    {
      this->i = i;
    }

    template<typename Scalar>
    VectorFormVol<Scalar>::VectorFormVol(unsigned int i, std::string area,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    VectorForm<Scalar>(i, area, ext, scaling_factor, u_ext_offset)
    {
    }

    template<typename Scalar>
    VectorFormVol<Scalar>::VectorFormVol(unsigned int i, Hermes::vector<std::string> areas,
      Hermes::vector<MeshFunction<Scalar>*> ext, double scaling_factor, int u_ext_offset) :
    VectorForm<Scalar>(i, areas, ext, scaling_factor, u_ext_offset)
    {
    }

    template<typename Scalar>
    Scalar VectorForm<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
      Geom<double> *e, ExtData<Scalar> *ext) const
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("VectorForm<Scalar>::value");
      return 0.0;
    }

    template<typename Scalar>
    Hermes::Ord VectorForm<Scalar>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("VectorForm<Scalar>::ord");
      return Hermes::Ord();
    }

    template<typename Scalar>
    VectorFormVol<Scalar>* VectorFormVol<Scalar>::clone()
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("VectorFormVol<Scalar>::clone()");
      return NULL;
    }

    template<typename Scalar>
    VectorFormSurf<Scalar>::VectorFormSurf(unsigned int i, std::string area,
      Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    VectorForm<Scalar>(i, area, ext, scaling_factor, u_ext_offset)
    {
    }
    template<typename Scalar>
    VectorFormSurf<Scalar>::VectorFormSurf(unsigned int i, Hermes::vector<std::string> areas,
      Hermes::vector<MeshFunction<Scalar>*> ext,
      double scaling_factor, int u_ext_offset) :
    VectorForm<Scalar>(i, areas, ext, scaling_factor, u_ext_offset)
    {
    }

    template<typename Scalar>
    VectorFormSurf<Scalar>* VectorFormSurf<Scalar>::clone()
    {
      throw Hermes::Exceptions::FunctionNotOverridenException("VectorFormSurf<Scalar>::clone()");
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
      seq++;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_matrix_form_surf(MatrixFormSurf<Scalar>* form)
    {
      if(form->i >= neq || form->j >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      mfsurf.push_back(form);
      seq++;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_vector_form(VectorFormVol<Scalar>* form)
    {
      if(form->i >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");
      form->set_weakform(this);
      vfvol.push_back(form);
      seq++;
    }

    template<typename Scalar>
    void WeakForm<Scalar>::add_vector_form_surf(VectorFormSurf<Scalar>* form)
    {
      if(form->i >= neq)
        throw Hermes::Exceptions::Exception("Invalid equation number.");

      form->set_weakform(this);
      vfsurf.push_back(form);
      seq++;
    }

    template<typename Scalar>
    Hermes::vector<MatrixFormVol<Scalar> *> WeakForm<Scalar>::get_mfvol()
    {
      return mfvol;
    }
    template<typename Scalar>
    Hermes::vector<MatrixFormSurf<Scalar> *> WeakForm<Scalar>::get_mfsurf()
    {
      return mfsurf;
    }
    template<typename Scalar>
      Hermes::vector<VectorFormVol<Scalar> *> WeakForm<Scalar>::get_vfvol()
    {
      return vfvol;
    }
    template<typename Scalar>
      Hermes::vector<VectorFormSurf<Scalar> *> WeakForm<Scalar>::get_vfsurf()
    {
      return vfsurf;
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
    template class HERMES_API MatrixFormVol<double>;
    template class HERMES_API MatrixFormVol<std::complex<double> >;
    template class HERMES_API MatrixFormSurf<double>;
    template class HERMES_API MatrixFormSurf<std::complex<double> >;
    template class HERMES_API VectorFormVol<double>;
    template class HERMES_API VectorFormVol<std::complex<double> >;
    template class HERMES_API VectorFormSurf<double>;
    template class HERMES_API VectorFormSurf<std::complex<double> >;
  }
}