// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "weakform_library/weakforms_neutronics.h"
#include "weakform_library/integrals_h1.h"

#include <algorithm>
#include <iomanip>
namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsNeutronics
    {
      namespace Monoenergetic
      {
        namespace Diffusion
        {
          template<typename Scalar>
          DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource( Hermes::vector<std::string> regions,
            Hermes::vector<double> D_map,
            Hermes::vector<double> Sigma_a_map,
            Hermes::vector<double> Q_map ) : WeakForm<Scalar>(1)
          {
            using namespace WeakFormsH1;

            for (unsigned int i = 0; i < regions.size(); i++)
            {
              /* Jacobian */
              // Diffusion.
              this->add_matrix_form(new DefaultJacobianDiffusion<Scalar>(0, 0, regions[i], new Hermes1DFunction<Scalar>(D_map[i]),
                HERMES_SYM));
              // Absorption.
              this->add_matrix_form(new DefaultMatrixFormVol<Scalar>(0, 0, regions[i], new Hermes2DFunction<Scalar>(Sigma_a_map[i]),
                HERMES_SYM));

              /* Residual */
              // Diffusion.
              this->add_vector_form(new DefaultResidualDiffusion<Scalar>(0, regions[i], new Hermes1DFunction<Scalar>(D_map[i])));
              // Absorption.
              this->add_vector_form(new DefaultResidualVol<Scalar>(0, regions[i], new Hermes2DFunction<Scalar>(Sigma_a_map[i])));
              // Sources.
              this->add_vector_form(new DefaultVectorFormVol<Scalar>(0, regions[i], new Hermes2DFunction<Scalar>(-Q_map[i])));
            }
          }
        }
      }

      namespace Multigroup
      {
        namespace MaterialProperties
        {
          namespace Common
          {
            void MaterialPropertyMaps::extend_to_multigroup(const MaterialPropertyMap0& mrsg_map,
              MaterialPropertyMap1 *mrmg_map)
            {
              if(G == 1)
                this->warn(W_MG_EXTENSION);

              MaterialPropertyMap0::const_iterator it;
              for (it = mrsg_map.begin(); it != mrsg_map.end(); ++it)
                (*mrmg_map)[it->first].assign(G, it->second);
            }

            void MaterialPropertyMaps::extend_to_multiregion(const rank1& srmg_array,
              MaterialPropertyMap1 *mrmg_map)
            {
              if(materials_list.empty())
                throw Hermes::Exceptions::Exception(E_MR_EXTENSION);

              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                (*mrmg_map)[*it] = srmg_array;
            }

            void MaterialPropertyMaps::extend_to_multiregion_multigroup(const rank0& srsg_value,
              MaterialPropertyMap1 *mrmg_map)
            {
              if(materials_list.empty())
                throw Hermes::Exceptions::Exception(E_MR_EXTENSION);
              if(G == 1)
                this->warn(W_MG_EXTENSION);

              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                (*mrmg_map)[*it].assign(G, srsg_value);
            }

            void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap1 *mrmg_map)
            {
              if(materials_list.empty())
                throw Hermes::Exceptions::Exception(E_MR_EXTENSION);

              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                (*mrmg_map)[*it].assign(G, c);
            }

            void MaterialPropertyMaps::validate()
            {
              using namespace ValidationFunctors;

              if(fission_multigroup_structure.empty())
                fission_multigroup_structure = bool1(G, true);

              if(chi.empty())
              {
                fill_with(0.0, &chi);
                MaterialPropertyMap1::iterator it = chi.begin();
                for ( ; it != chi.end(); ++it)
                  it->second[0] = 1.0;
                fission_multigroup_structure = bool1(G, false);
                fission_multigroup_structure[0] = true;
              }

              if(nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
                nu = NDArrayMapOp::divide<rank1>(nuSigma_f, Sigma_f);
              else if(nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
                nuSigma_f = NDArrayMapOp::multiply<rank1>(nu, Sigma_f);
              else if(Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
                Sigma_f = NDArrayMapOp::divide<rank1>(nuSigma_f, nu);
              else if(!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
              {
                MaterialPropertyMap1 diff = NDArrayMapOp::subtract<rank1>(nuSigma_f,
                  NDArrayMapOp::multiply<rank1>(nu, Sigma_f) );
                std::for_each(diff.begin(), diff.end(), ensure_trivial());
              }
              else
              {
                this->warn(W_NO_FISSION);
                fill_with(0.0, &nu);
                fill_with(0.0, &chi);
                fill_with(0.0, &Sigma_f);
              }

              if((nu.size() != Sigma_f.size()) || (nu.size() != chi.size()))
                throw Hermes::Exceptions::Exception(E_NONMATCHING_PROPERTIES);

              if(Sigma_f.size() > 0)
              {
                std::for_each(nu.begin(), nu.end(), ensure_size(G));
                std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size(G));
                std::for_each(chi.begin(), chi.end(), ensure_size(G));
              }

              if(Sigma_a.size() > 0)
              {
                // Warn if \Sigma_a < \Sigma_f for any region (this indicates an unphysical situation, since
                // by definition \Sigma_a = \Sigma_f + \Sigma_c + \Sigma_{n, p} + other possible reactions
                // leading to neutron removal).
                MaterialPropertyMap1::const_iterator ita = Sigma_a.begin();
                MaterialPropertyMap1::const_iterator itf = Sigma_f.begin();
                for ( ; ita != Sigma_a.end(); ++ita, ++itf)
                {
                  rank1::const_iterator a = ita->second.begin();
                  rank1::const_iterator f = itf->second.begin();

                  for ( ; a != ita->second.end(); ++a, ++f)
                    if(*a < *f)
                      this->warn(W_SA_LT_SF);
                }
              }
            }

            const rank1& MaterialPropertyMaps::get_Sigma_f(std::string material) const
            {
              if(material == "-999") return this->Sigma_f.begin()->second;

              // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
              // which operator[] is undefined.
              MaterialPropertyMap1::const_iterator data = this->Sigma_f.find(material);
              if(data != this->Sigma_f.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }
            const rank1& MaterialPropertyMaps::get_nu(std::string material) const
            {
              if(material == "-999") return this->nu.begin()->second;

              MaterialPropertyMap1::const_iterator data = this->nu.find(material);
              if(data != this->nu.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }
            const rank1& MaterialPropertyMaps::get_chi(std::string material) const
            {
              if(material == "-999") return this->chi.begin()->second;

              MaterialPropertyMap1::const_iterator data = this->chi.find(material);
              if(data != this->chi.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }

            std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
            {
              using namespace std;

              os << std::endl;
              os << setw(12) << "target group" << setw(10) << "chi" << setw(10) << "nu";
              os << setw(10) << "Sigma_f" << std::endl;

              MaterialPropertyMap1::const_iterator data_elem = matprop.chi.begin();
              for ( ; data_elem != matprop.chi.end(); ++data_elem)
              {
                string mat = data_elem->first;

                os << setw(80) << setfill('-') << ' ' << std::endl << setfill(' ');
                os << setw(40) << mat << std::endl;
                os << setw(80) << setfill('-') << ' ' << std::endl << setfill(' ');
                for (unsigned int gto = 0; gto < matprop.G; gto++)
                {
                  os << setw(6) << gto << setw(6) << ' ';
                  os << setw(10) << matprop.get_chi(mat)[gto];
                  os << setw(10) << matprop.get_nu(mat)[gto];
                  os << setw(10) << matprop.get_Sigma_f(mat)[gto];

                  os << std::endl;
                }
              }

              os << std::endl;
              return os;
            }
          }

          namespace Diffusion
          {
            MaterialPropertyMap1 MaterialPropertyMaps::extract_map2_diagonals(const MaterialPropertyMap2& map2)
            {
              MaterialPropertyMap1 diags;

              MaterialPropertyMap2::const_iterator map2_it = map2.begin();
              for ( ; map2_it != map2.end(); ++map2_it)
              {
                diags[map2_it->first].reserve(G);
                for (unsigned int g = 0; g < G; g++)
                  diags[map2_it->first].push_back(map2_it->second[g][g]);
              }

              return diags;
            }

            MaterialPropertyMap1 MaterialPropertyMaps::sum_map2_columns(const MaterialPropertyMap2& map2)
            {
              MaterialPropertyMap1 summed;

              MaterialPropertyMap2::const_iterator map2_it = map2.begin();
              for ( ; map2_it != map2.end(); ++map2_it)
              {
                summed[map2_it->first].reserve(G);
                for (unsigned int gfrom = 0; gfrom < G; gfrom++)
                {
                  double sum = 0.0;

                  for (unsigned int gto = 0; gto < G; gto++)
                    sum += map2_it->second[gto][gfrom];

                  summed[map2_it->first].push_back(sum);
                }
              }

              return summed;
            }

            MaterialPropertyMap2 MaterialPropertyMaps::create_map2_by_diagonals(const MaterialPropertyMap1& diags)
            {
              MaterialPropertyMap2 map2;

              MaterialPropertyMap1::const_iterator diags_it = diags.begin();
              for ( ; diags_it != diags.end(); ++diags_it)
              {
                map2[diags_it->first].resize(G, rank1(G, 0.0));

                for (unsigned int g = 0; g < G; g++)
                  map2[diags_it->first][g][g] = diags_it->second[g];
              }

              return map2;
            }

            void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap2 *mrmg_map)
            {
              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                (*mrmg_map)[*it].assign(G, rank1(G, c));
            }

            void MaterialPropertyMaps::validate()
            {
              Common::MaterialPropertyMaps::validate();

              bool D_given = !D.empty();
              bool Sigma_r_given = !Sigma_r.empty();
              bool Sigma_s_given = !Sigma_s.empty();
              bool Sigma_t_given = !Sigma_t.empty();
              bool Sigma_a_given = !Sigma_a.empty();
              bool Sigma_f_given = !Sigma_f.empty();
              bool src_given = !src.empty();

              if(!Sigma_r_given)
              {
                // If Sigma_r is not given, we can calculate it from Sigma_t and Sigma_s.

                if(Sigma_t_given)
                {
                  if(!Sigma_s_given)
                  {
                    if(Sigma_a_given)
                    {
                      // If Sigma_s is not given, but Sigma_a is, we can calculate Sigma_s from Sigma_t and Sigma_a.
                      Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_a));
                    }
                    else
                    {
                      // If only Sigma_t is given, we assume that all reaction terms are included in Sigma_t; all
                      // other x-sections will be set to zero.
                      this->warn(W_NO_SCATTERING);
                      fill_with(0.0, &Sigma_s);
                    }

                    Sigma_s_given = true;
                  }
                }
                else
                {
                  // If Sigma_t is not given, but Sigma_a and Sigma_s are, we can obtain Sigma_t from the latter two.

                  if(!Sigma_s_given)
                  {
                    this->warn(W_NO_SCATTERING);
                    fill_with(0.0, &Sigma_s);
                    Sigma_s_given = true;
                  }

                  if(Sigma_a_given)
                    Sigma_t = Common::NDArrayMapOp::add<rank1>(Sigma_a, sum_map2_columns(Sigma_s));
                  else
                  {
                    // If neither Sigma_r, Sigma_t, Sigma_a are given, we may have a purely fissioning system.
                    if(Sigma_f_given)
                      Sigma_t = Sigma_f;
                    else
                      throw Hermes::Exceptions::Exception(E_INSUFFICIENT_DATA);
                  }

                  Sigma_t_given = true;
                }

                Sigma_r = Common::NDArrayMapOp::subtract<rank1>(Sigma_t, extract_map2_diagonals(Sigma_s));
                Sigma_r_given = true;
              }

              // Now, we surely have Sigma_r ...

              if(scattering_multigroup_structure.empty())
                scattering_multigroup_structure = bool2(G, Hermes::vector<bool>(G, true));

              if(!Sigma_s_given)
              {
                // If Sigma_s is not given, but Sigma_t is, we can obtain the former from the latter and from Sigma_r.
                // Note that the execution will come here only if the user entered Sigma_r himself - otherwise, Sigma_s
                // has been already set in the previous test case.

                if(Sigma_t_given)
                {
                  Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_r));

                  scattering_multigroup_structure = bool2(G, Hermes::vector<bool>(G, false));
                  for (unsigned int gto = 0; gto < G; gto++)
                    for (unsigned int gfrom = 0; gfrom < G; gfrom++)
                      if(gto == gfrom)
                        scattering_multigroup_structure[gto][gfrom] = true;
                }
                else
                {
                  this->warn(W_NO_SCATTERING);
                  fill_with(0.0, &Sigma_s);
                  scattering_multigroup_structure = bool2(G, Hermes::vector<bool>(G, false));
                }

                Sigma_s_given = true;
              }

              // Now, we surely have Sigma_s and Sigma_r, one parameter to go ...

              if(!D_given)
              {
                MaterialPropertyMap1::const_iterator Sr_elem = Sigma_r.begin();
                for ( ; Sr_elem != Sigma_r.end(); ++Sr_elem)
                  for (unsigned int g = 0; g < G; g++)
                    D[Sr_elem->first][g] = 1./(3.*Sr_elem->second[g]);

                D_given = true;
              }

              if((D.size() != Sigma_r.size()) || (D.size() != Sigma_s.size()) || (src_given && D.size() != src.size()))
                throw Hermes::Exceptions::Exception(E_NONMATCHING_PROPERTIES);

              using ValidationFunctors::ensure_size;
              std::for_each(Sigma_s.begin(), Sigma_s.end(), ensure_size(G, G));
              std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size(G));
              std::for_each(src.begin(), src.end(), ensure_size(G));
              std::for_each(D.begin(), D.end(), ensure_size(G));
            }

            const rank2& MaterialPropertyMaps::get_Sigma_s(std::string material) const
            {
              if(material == "-999") return this->Sigma_s.begin()->second;

              // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
              // which operator[] is undefined.
              MaterialPropertyMap2::const_iterator data = this->Sigma_s.find(material);
              if(data != this->Sigma_s.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank2()); // To avoid MSVC problems; execution should never come to this point.
              }
            }
            const rank1& MaterialPropertyMaps::get_Sigma_r(std::string material) const
            {
              if(material == "-999") return this->Sigma_r.begin()->second;

              MaterialPropertyMap1::const_iterator data = this->Sigma_r.find(material);
              if(data != this->Sigma_r.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }
            const rank1& MaterialPropertyMaps::get_D(std::string material) const
            {
              if(material == "-999") return this->D.begin()->second;

              MaterialPropertyMap1::const_iterator data = this->D.find(material);
              if(data != this->D.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }
            const rank1& MaterialPropertyMaps::get_src(std::string material) const
            {
              if(material == "-999") return this->src.begin()->second;

              MaterialPropertyMap1::const_iterator data = this->src.find(material);
              if(data != this->src.end())
                return data->second;
              else
              {
                throw Hermes::Exceptions::Exception(E_INVALID_MARKER);
                return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
              }
            }

            HERMES_API std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
            {
              using namespace std;

              os << static_cast<const Common::MaterialPropertyMaps&>(matprop) << std::endl;

              os << setw(12) << "target group" << setw(10) << "D" << setw(10) << "Sigma_r";
              os << setw(10) << "ext. src" << setw(22) << "Sigma_s" << std::endl;

              MaterialPropertyMap1::const_iterator data_elem = matprop.Sigma_r.begin();
              for ( ; data_elem != matprop.Sigma_r.end(); ++data_elem)
              {
                string mat = data_elem->first;

                os << setw(80) << setfill('-') << ' ' << std::endl << setfill(' ');
                os << setw(40) << mat << std::endl;
                os << setw(80) << setfill('-') << ' ' << std::endl << setfill(' ');
                for (unsigned int gto = 0; gto < matprop.G; gto++)
                {
                  os << setw(6) << gto << setw(6) << ' ';
                  os << setw(10) << matprop.get_D(mat)[gto];
                  os << setw(10) << matprop.get_Sigma_r(mat)[gto];
                  os << setw(10);
                  if(matprop.src.empty())
                    os << "N/A";
                  else
                    os << matprop.get_src(mat)[gto];

                  for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
                    os << setw(8) << matprop.get_Sigma_s(mat)[gto][gfrom];

                  os << std::endl;
                }
              }

              return os << std::endl;
            }
          }
        }

        namespace ElementaryForms
        {
          namespace Diffusion
          {
            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar VacuumBoundaryCondition::Jacobian<ScalarClass>::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              Scalar result;

              if(geom_type == HERMES_PLANAR)
                result = 0.5 * int_u_v<Real, Scalar>(n, wt, u, v);
              else if(geom_type == HERMES_AXISYM_X)
                result = 0.5 * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
              else
                result = 0.5 * int_x_u_v<Real, Scalar>(n, wt, u, v, e);

              return result;
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar VacuumBoundaryCondition::Residual<ScalarClass>::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              Scalar result;

              if(geom_type == HERMES_PLANAR)
                result = 0.5 * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
              else if(geom_type == HERMES_AXISYM_X)
                result = 0.5 * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
              else
                result = 0.5 * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);

              return result;
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar DiffusionReaction::Jacobian<ScalarClass>::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              Scalar result;

              std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
              rank1 D_elem = matprop.get_D(mat);
              rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);

              if(geom_type == HERMES_PLANAR)
              {
                result = D_elem[g] * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
                  Sigma_r_elem[g] * int_u_v<Real, Scalar>(n, wt, u, v);
              }
              else
              {
                if(geom_type == HERMES_AXISYM_X)
                {
                  result = D_elem[g] * int_y_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
                    Sigma_r_elem[g] * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
                }
                else
                {
                  result = D_elem[g] * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
                    Sigma_r_elem[g] * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
                }
              }
              return result;
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar DiffusionReaction::Residual<ScalarClass>::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              Scalar result;

                            std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
              rank1 D_elem = matprop.get_D(mat);
              rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);

              if(geom_type == HERMES_PLANAR)
              {
                result = D_elem[g] * int_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v) +
                  Sigma_r_elem[g] * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
              }
              else
              {
                if(geom_type == HERMES_AXISYM_X)
                {
                  result = D_elem[g] * int_y_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) +
                    Sigma_r_elem[g] * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
                }
                else
                {
                  result = D_elem[g] * int_x_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) +
                    Sigma_r_elem[g] * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
                }
              }
              return result;
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar FissionYield::Jacobian<ScalarClass>::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              if(!matprop.get_fission_multigroup_structure()[gto])
                return Scalar(0.0);

              Scalar result = Scalar(0);
              if(geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
              else
              {
                if(geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
                else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
              }

              std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
              rank1 nu_elem = matprop.get_nu(mat);
              rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
              rank1 chi_elem = matprop.get_chi(mat);

              return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar FissionYield::OuterIterationForm<ScalarClass>::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
              if(!matprop.get_fission_multigroup_structure()[g])
                return Scalar(0.0);

              std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
              rank1 nu_elem = matprop.get_nu(mat);
              rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
              rank1 chi_elem = matprop.get_chi(mat);

              Scalar result = Scalar(0);
              for (int i = 0; i < n; i++)
              {
                Scalar local_res = Scalar(0);
                for (int gfrom = 0; gfrom < this->wf->get_ext().size(); gfrom++)
                  local_res += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext[gfrom]->val[i];

                local_res = local_res * wt[i] * v->val[i];

                if(geom_type == HERMES_AXISYM_X)
                  local_res = local_res * e->y[i];
                else if(geom_type == HERMES_AXISYM_Y)
                  local_res = local_res * e->x[i];

                result += local_res;
              }

              return result * chi_elem[g] / keff;
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar FissionYield::Residual<ScalarClass>::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext ) const
            {
              if(!matprop.get_fission_multigroup_structure()[gto])
                return Scalar(0.0);

              Scalar result = Scalar(0);
              if(geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
              else
              {
                if(geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
              }

              std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
              rank1 nu_elem = matprop.get_nu(mat);
              rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
              rank1 chi_elem = matprop.get_chi(mat);

              return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar Scattering::Jacobian<ScalarClass>::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext  ) const
            {
              Scalar result = Scalar(0);
              if(geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
              else
              {
                if(geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
                else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
              }

              return result * matprop.get_Sigma_s(this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker)[gto][gfrom];
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar Scattering::Residual<ScalarClass>::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext ) const
            {
              Scalar result = Scalar(0);
              if(geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
              else
              {
                if(geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
              }

              return result * matprop.get_Sigma_s(this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker)[gto][gfrom];
            }

            template<typename ScalarClass>
            template<typename Real, typename Scalar>
            Scalar ExternalSources::LinearForm<ScalarClass>::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
              Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const
            {
                            std::string mat = this->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;

              if(geom_type == HERMES_PLANAR)
                return matprop.get_src(mat)[g] * int_v<Real>(n, wt, v);
              else
              {
                if(geom_type == HERMES_AXISYM_X)
                  return matprop.get_src(mat)[g] * int_y_v<Real>(n, wt, v, e);
                else
                  return matprop.get_src(mat)[g] * int_x_v<Real>(n, wt, v, e);
              }
            }
          }
        }

        namespace CompleteWeakForms
        {
          namespace Diffusion
          {
            template<typename Scalar>
            void DefaultWeakFormFixedSource<Scalar>::lhs_init(unsigned int G, const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              GeomType geom_type)
            {
              bool2 Ss_nnz = matprop.get_scattering_multigroup_structure();
              bool1 chi_nnz = matprop.get_fission_multigroup_structure();

              for (unsigned int gto = 0; gto < G; gto++)
              {
                this->add_matrix_form(new DiffusionReaction::Jacobian<Scalar>(gto, matprop, mesh, geom_type));
                this->add_vector_form(new DiffusionReaction::Residual<Scalar>(gto, matprop, mesh, geom_type));

                for (unsigned int gfrom = 0; gfrom < G; gfrom++)
                {
                  if(Ss_nnz[gto][gfrom])
                  {
                    this->add_matrix_form(new Scattering::Jacobian<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                    this->add_vector_form(new Scattering::Residual<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                  }

                  if(chi_nnz[gto])
                  {
                    this->add_matrix_form(new FissionYield::Jacobian<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                    this->add_vector_form(new FissionYield::Residual<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                  }
                }
              }
            }

            template<typename Scalar>
            DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              GeomType geom_type) : WeakForm<Scalar>(matprop.get_G())
            {
              lhs_init(matprop.get_G(), matprop, mesh, geom_type);
              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
                this->add_vector_form(new ExternalSources::LinearForm<Scalar>(gto, matprop, mesh, geom_type));
            }

            template<typename Scalar>
            DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              Hermes2DFunction<Scalar>*f_src, std::string src_area,
              GeomType geom_type  ) : WeakForm<Scalar>(matprop.get_G())
            {
              lhs_init(matprop.get_G(), matprop, mesh, geom_type);
              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
                this->add_vector_form(new WeakFormsH1::DefaultVectorFormVol<Scalar>(gto, src_area, f_src, geom_type));
            }

            template<typename Scalar>
            DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              Hermes2DFunction<Scalar>*f_src,
              Hermes::vector<std::string> src_areas,
              GeomType geom_type  ) : WeakForm<Scalar>(matprop.get_G())
            {
              lhs_init(matprop.get_G(), matprop, mesh, geom_type);
              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
                this->add_vector_form(new WeakFormsH1::DefaultVectorFormVol<Scalar>(gto, src_areas, f_src, geom_type));
            }

            template<typename Scalar>
            DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              const Hermes::vector<Hermes2DFunction<Scalar>*>& f_src,
              std::string src_area,
              GeomType geom_type ) : WeakForm<Scalar>(matprop.get_G())
            {
              if(f_src.size() != matprop.get_G())
                throw Hermes::Exceptions::Exception(E_INVALID_SIZE);

              lhs_init(matprop.get_G(), matprop, mesh, geom_type);
              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
                this->add_vector_form(new WeakFormsH1::DefaultVectorFormVol<Scalar>(gto, src_area, f_src[gto], geom_type));
            }

            template<typename Scalar>
            DefaultWeakFormFixedSource<Scalar>::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              const Hermes::vector<Hermes2DFunction<Scalar>*>& f_src,
              Hermes::vector<std::string> src_areas,
              GeomType geom_type ) : WeakForm<Scalar>(matprop.get_G())
            {
              if(f_src.size() != matprop.get_G())
                throw Hermes::Exceptions::Exception(E_INVALID_SIZE);

              lhs_init(matprop.get_G(), matprop, mesh, geom_type);
              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
                this->add_vector_form(new WeakFormsH1::DefaultVectorFormVol<Scalar>(gto, src_areas, f_src[gto], geom_type));
            }

            template<typename Scalar>
            DefaultWeakFormSourceIteration<Scalar>::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop, MeshSharedPtr mesh,
              Hermes::vector<MeshFunctionSharedPtr<Scalar> >& iterates,
              double initial_keff_guess,
              GeomType geom_type ) : WeakForm<Scalar>(matprop.get_G())
            {
              bool2 Ss_nnz = matprop.get_scattering_multigroup_structure();

              for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
              {
                this->add_matrix_form(new DiffusionReaction::Jacobian<Scalar>(gto, matprop, mesh, geom_type));
                this->add_vector_form(new DiffusionReaction::Residual<Scalar>(gto, matprop, mesh, geom_type));

                for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
                {
                  if(Ss_nnz[gto][gfrom])
                  {
                    this->add_matrix_form(new Scattering::Jacobian<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                    this->add_vector_form(new Scattering::Residual<Scalar>(gto, gfrom, matprop, mesh, geom_type));
                  }
                }

                FissionYield::OuterIterationForm<Scalar>* keff_iteration_form =
                  new FissionYield::OuterIterationForm<Scalar>( gto, matprop, mesh, iterates, initial_keff_guess, geom_type );
                keff_iteration_forms.push_back(keff_iteration_form);
                this->add_vector_form(keff_iteration_form);
              }
            }

            template<typename Scalar>
            void DefaultWeakFormSourceIteration<Scalar>::update_keff(double new_keff)
            {
              /* Somehow does not work with templates. A bug / typo from me?
              Hermes::vector<FissionYield::OuterIterationForm<Scalar> *>::iterator it = keff_iteration_forms.begin();
              for ( ; it != keff_iteration_forms.end(); ++it)
              (*it)->update_keff(new_keff);
              */
              for(int i = 0; i < keff_iteration_forms.size(); i++)
                keff_iteration_forms[i]->update_keff(new_keff);
            }
          }
        }

        namespace SupportClasses
        {
          void SourceFilter::filter_fn(int n, Hermes::vector<double*> values, double* result)
          {
            for (int i = 0; i < n; i++)
            {
              result[i] = 0;
              for (unsigned int j = 0; j < values.size(); j++)
                result[i] += nu[j] * Sigma_f[j] * values.at(j)[i];
            }
          }
        }
      }
      namespace Monoenergetic
      {
        namespace Diffusion
        {
          template class HERMES_API DefaultWeakFormFixedSource<double>;
          template class HERMES_API DefaultWeakFormFixedSource<std::complex<double> >;
        }
      }
      namespace Multigroup
      {
        namespace ElementaryForms
        {
          namespace Diffusion
          {
            template class HERMES_API FissionYield::Jacobian<double>;
            template class HERMES_API FissionYield::Jacobian<std::complex<double> >;
            template class HERMES_API FissionYield::OuterIterationForm<double>;
            template class HERMES_API FissionYield::OuterIterationForm<std::complex<double> >;
            template class HERMES_API FissionYield::Residual<double>;
            template class HERMES_API FissionYield::Residual<std::complex<double> >;

            template class HERMES_API VacuumBoundaryCondition::Jacobian<double>;
            template class HERMES_API VacuumBoundaryCondition::Jacobian<std::complex<double> >;
            template class HERMES_API VacuumBoundaryCondition::Residual<double>;
            template class HERMES_API VacuumBoundaryCondition::Residual<std::complex<double> >;

            template double VacuumBoundaryCondition::Jacobian<double>::matrix_form<double, double>(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;

            template double VacuumBoundaryCondition::Residual<double>::vector_form<double, double>(int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;

            template double DiffusionReaction::Jacobian<double>::matrix_form<double, double>(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;

            template double DiffusionReaction::Residual<double>::vector_form<double, double>(int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;

            template double FissionYield::Jacobian<double>::matrix_form<double, double>( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;

            template double FissionYield::OuterIterationForm<double>::vector_form<double, double>( int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext ) const;

            template double FissionYield::Residual<double>::vector_form<double, double>( int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext ) const;

            template double Scattering::Jacobian<double>::matrix_form<double, double>( int n, double *wt, Func<double> *u_ext[], Func<double> *u,
              Func<double> *v, Geom<double> *e, Func<double> **ext  ) const;

            template double Scattering::Residual<double>::vector_form<double, double>( int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext ) const;

            template double ExternalSources::LinearForm<double>::vector_form<double, double>(int n, double *wt, Func<double> *u_ext[],
              Func<double> *v, Geom<double> *e, Func<double> **ext) const;
          }
        }

        namespace CompleteWeakForms
        {
          namespace Diffusion
          {
            template class HERMES_API DefaultWeakFormFixedSource<double>;
            template class HERMES_API DefaultWeakFormFixedSource<std::complex<double> >;

            template class HERMES_API DefaultWeakFormSourceIteration<double>;
            template class HERMES_API DefaultWeakFormSourceIteration<std::complex<double> >;
          }
        }
      }
    }
  }
}
