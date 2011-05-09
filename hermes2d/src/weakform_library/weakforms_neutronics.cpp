#include "../hermes2d.h"

#include <algorithm>
#include <iomanip>

namespace WeakFormsNeutronics
{
  namespace Monoenergetic
  {    
    namespace Diffusion 
    {
      DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( Hermes::vector<std::string> regions, 
                                                              Hermes::vector<double> D_map, 
                                                              Hermes::vector<double> Sigma_a_map, 
                                                              Hermes::vector<double> Q_map ) : WeakForm(1) 
      {
        using namespace WeakFormsH1;
        
        for (unsigned int i = 0; i < regions.size(); i++)
        {
          /* Jacobian */
          // Diffusion.
          add_matrix_form(new DefaultJacobianDiffusion(0, 0, regions[i], D_map[i], 
                                                       HERMES_DEFAULT_SPLINE, HERMES_SYM));
          // Absorption.
          add_matrix_form(new DefaultMatrixFormVol(0, 0, regions[i], Sigma_a_map[i], 
                                                   HERMES_DEFAULT_FUNCTION, HERMES_SYM));
          
          /* Residual */
          // Diffusion.
          add_vector_form(new DefaultResidualDiffusion(0, regions[i], D_map[i]));
          // Absorption.
          add_vector_form(new DefaultResidualVol(0, regions[i], Sigma_a_map[i]));
          // Sources.
          add_vector_form(new DefaultVectorFormVol(0, regions[i], -Q_map[i]));
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
          if (G == 1)
            warning(W_MG_EXTENSION);
          
          MaterialPropertyMap0::const_iterator it;
          for (it = mrsg_map.begin(); it != mrsg_map.end(); ++it)
            (*mrmg_map)[it->first].assign(G, it->second);
          
        }
        
        void MaterialPropertyMaps::extend_to_multiregion(const rank1& srmg_array, 
                                                         MaterialPropertyMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            (*mrmg_map)[*it] = srmg_array;
        }
        
        void MaterialPropertyMaps::extend_to_multiregion_multigroup(const rank0& srsg_value, 
                                                                    MaterialPropertyMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          if (G == 1)
            warning(W_MG_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            (*mrmg_map)[*it].assign(G, srsg_value);
        }
        
        void MaterialPropertyMaps::fill_with(double c, MaterialPropertyMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            (*mrmg_map)[*it].assign(G, c);
        }
        
        void MaterialPropertyMaps::validate()
        {       
          using namespace ValidationFunctors;
          
          if (fission_multigroup_structure.empty())
            fission_multigroup_structure = bool1(G, true);
          
          if (chi.empty())
          {
            fill_with(0.0, &chi);
            MaterialPropertyMap1::iterator it = chi.begin();
            for ( ; it != chi.end(); ++it)
              it->second[0] = 1.0;
            fission_multigroup_structure = bool1(G, false);
            fission_multigroup_structure[0] = true;
          }
          
          if (nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
            nu = NDArrayMapOp::divide<rank1>(nuSigma_f, Sigma_f);
          else if (nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
            nuSigma_f = NDArrayMapOp::multiply<rank1>(nu, Sigma_f);
          else if (Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
            Sigma_f = NDArrayMapOp::divide<rank1>(nuSigma_f, nu);
          else if (!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
          {
            MaterialPropertyMap1 diff = NDArrayMapOp::subtract<rank1>(nuSigma_f, 
                                                                      NDArrayMapOp::multiply<rank1>(nu, Sigma_f) );
            std::for_each(diff.begin(), diff.end(), ensure_trivial());
          }
          else
          {
            warning(W_NO_FISSION);
            fill_with(0.0, &nu);
            fill_with(0.0, &chi);
            fill_with(0.0, &Sigma_f);
          }
          
          if ((nu.size() != Sigma_f.size()) || (nu.size() != chi.size()))
            error(E_NONMATCHING_PROPERTIES);
          
          if (Sigma_f.size() > 0)
          {
            std::for_each(nu.begin(), nu.end(), ensure_size(G));
            std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size(G));
            std::for_each(chi.begin(), chi.end(), ensure_size(G));
          }
          
          if (Sigma_a.size() > 0)
          {
            // Warn if \Sigma_a < \Sigma_f for any region (this indicates an unphysical situation, since
            // by definition \Sigma_a = \Sigma_f + \Sigma_c + \Sigma_{n,p} + other possible reactions
            // leading to neutron removal).
            MaterialPropertyMap1::const_iterator ita = Sigma_a.begin();
            MaterialPropertyMap1::const_iterator itf = Sigma_f.begin();
            for ( ; ita != Sigma_a.end(); ++ita, ++itf)
            {
              rank1::const_iterator a = ita->second.begin();
              rank1::const_iterator f = itf->second.begin();
              
              for ( ; a != ita->second.end(); ++a,++f)
                if (*a < *f)
                  warning(W_SA_LT_SF);
            }
          }
        }
        
        const rank1& MaterialPropertyMaps::get_Sigma_f(std::string material) const
        {
          // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
          // which operator[] is undefined.
          MaterialPropertyMap1::const_iterator data = this->Sigma_f.find(material);
          if (data != this->Sigma_f.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        const rank1& MaterialPropertyMaps::get_nu(std::string material) const
        {
          MaterialPropertyMap1::const_iterator data = this->nu.find(material);
          if (data != this->nu.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        const rank1& MaterialPropertyMaps::get_chi(std::string material) const
        {
          MaterialPropertyMap1::const_iterator data = this->chi.find(material);
          if (data != this->chi.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        
        std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
        {
          using namespace std;
          
          os << endl;
          os << setw(12) << "target group" << setw(10) << "chi" << setw(10) << "nu";
          os << setw(10) << "Sigma_f" << endl; 
          
          MaterialPropertyMap1::const_iterator data_elem = matprop.chi.begin();
          for ( ; data_elem != matprop.chi.end(); ++data_elem)
          {
            string mat = data_elem->first;
            
            os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
            os << setw(40) << mat << endl;
            os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
            for (unsigned int gto = 0; gto < matprop.G; gto++)
            {
              os << setw(6) << gto << setw(6) << ' ';
              os << setw(10) << matprop.get_chi(mat)[gto];
              os << setw(10) << matprop.get_nu(mat)[gto];
              os << setw(10) << matprop.get_Sigma_f(mat)[gto];
              
              os << endl;
            }
          }
          
          os << endl;
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
          
          if (!Sigma_r_given)
          {
            // If Sigma_r is not given, we can calculate it from Sigma_t and Sigma_s.
            
            if (Sigma_t_given)
            {
              if (!Sigma_s_given)
              {
                if (Sigma_a_given)
                {
                  // If Sigma_s is not given, but Sigma_a is, we can calculate Sigma_s from Sigma_t and Sigma_a.
                  Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_a));
                }
                else 
                {
                  // If only Sigma_t is given, we assume that all reaction terms are included in Sigma_t; all
                  // other x-sections will be set to zero.
                  warning(W_NO_SCATTERING);
                  fill_with(0.0, &Sigma_s);
                }
                
                Sigma_s_given = true;
              }
            }
            else
            {
              // If Sigma_t is not given, but Sigma_a and Sigma_s are, we can obtain Sigma_t from the latter two.
              
              if (!Sigma_s_given)
              {
                warning(W_NO_SCATTERING);
                fill_with(0.0, &Sigma_s);
                Sigma_s_given = true;
              }
              
              if (Sigma_a_given)
                Sigma_t = Common::NDArrayMapOp::add<rank1>(Sigma_a, sum_map2_columns(Sigma_s));
              else 
              {
                // If neither Sigma_r, Sigma_t, Sigma_a are given, we may have a purely fissioning system.
                if (Sigma_f_given)
                  Sigma_t = Sigma_f;
                else
                  error(E_INSUFFICIENT_DATA);
              }
              
              Sigma_t_given = true;
            }
            
            Sigma_r = Common::NDArrayMapOp::subtract<rank1>(Sigma_t, extract_map2_diagonals(Sigma_s));
            Sigma_r_given = true;
          }
          
          // Now, we surely have Sigma_r ...
          
          if (scattering_multigroup_structure.empty())
            scattering_multigroup_structure = bool2(G, std::vector<bool>(G, true));
          
          if (!Sigma_s_given)
          {
            // If Sigma_s is not given, but Sigma_t is, we can obtain the former from the latter and from Sigma_r.
            // Note that the execution will come here only if the user entered Sigma_r himself - otherwise, Sigma_s
            // has been already set in the previous test case.
            
            if (Sigma_t_given)
            {
              Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_r));
              
              scattering_multigroup_structure = bool2(G, std::vector<bool>(G, false));
              for (unsigned int gto = 0; gto < G; gto++)
                for (unsigned int gfrom = 0; gfrom < G; gfrom++)
                  if (gto == gfrom) 
                    scattering_multigroup_structure[gto][gfrom] = true;
            }
            else
            {
              warning(W_NO_SCATTERING);
              fill_with(0.0, &Sigma_s);
              scattering_multigroup_structure = bool2(G, std::vector<bool>(G, false));
            }
            
            Sigma_s_given = true;
          }
          
          // Now, we surely have Sigma_s and Sigma_r, one parameter to go ...
          
          if (!D_given)
          {
            MaterialPropertyMap1::const_iterator Sr_elem = Sigma_r.begin();
            for ( ; Sr_elem != Sigma_r.end(); ++Sr_elem)
              for (unsigned int g = 0; g < G; g++)
                D[Sr_elem->first][g] = 1./(3.*Sr_elem->second[g]);
              
              D_given = true;
          }
          
          if ((D.size() != Sigma_r.size()) || (D.size() != Sigma_s.size()) || (src_given && D.size() != src.size()))
            error(E_NONMATCHING_PROPERTIES);
          
          using ValidationFunctors::ensure_size;
          std::for_each(Sigma_s.begin(), Sigma_s.end(), ensure_size(G,G));
          std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size(G));
          std::for_each(src.begin(), src.end(), ensure_size(G));
          std::for_each(D.begin(), D.end(), ensure_size(G));
        }
        
        const rank2& MaterialPropertyMaps::get_Sigma_s(std::string material) const
        {
          // Note that prop[e->elem_marker] cannot be used since 'prop' is a constant std::map for
          // which operator[] is undefined.
          MaterialPropertyMap2::const_iterator data = this->Sigma_s.find(material);
          if (data != this->Sigma_s.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank2()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        const rank1& MaterialPropertyMaps::get_Sigma_r(std::string material) const
        {
          MaterialPropertyMap1::const_iterator data = this->Sigma_r.find(material);
          if (data != this->Sigma_r.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        const rank1& MaterialPropertyMaps::get_D(std::string material) const
        {
          MaterialPropertyMap1::const_iterator data = this->D.find(material);
          if (data != this->D.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        const rank1& MaterialPropertyMaps::get_src(std::string material) const
        {
          MaterialPropertyMap1::const_iterator data = this->src.find(material);
          if (data != this->src.end())
            return data->second;
          else
          {
            error(E_INVALID_MARKER);
            return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
          }
        }
        
        std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop)
        {
          using namespace std;
          
          os << static_cast<const Common::MaterialPropertyMaps&>(matprop) << endl;
          
          os << setw(12) << "target group" << setw(10) << "D" << setw(10) << "Sigma_r";
          os << setw(10) << "ext. src" << setw(22) << "Sigma_s" << endl; 
          
          MaterialPropertyMap1::const_iterator data_elem = matprop.Sigma_r.begin();
          for ( ; data_elem != matprop.Sigma_r.end(); ++data_elem)
          {
            string mat = data_elem->first;
            
            os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
            os << setw(40) << mat << endl;
            os << setw(80) << setfill('-') << ' ' << endl << setfill(' ');
            for (unsigned int gto = 0; gto < matprop.G; gto++)
            {
              os << setw(6) << gto << setw(6) << ' ';
              os << setw(10) << matprop.get_D(mat)[gto];
              os << setw(10) << matprop.get_Sigma_r(mat)[gto];
              os << setw(10);
              if (matprop.src.empty())
                os << "N/A";
              else
                os << matprop.get_src(mat)[gto];
              
              for (unsigned int gfrom = 0; gfrom < matprop.G; gfrom++)
                os << setw(8) << matprop.get_Sigma_s(mat)[gto][gfrom];
              
              os << endl;
            }
          }
          
          return os << endl;
        }
      }
    }
    
    namespace ElementaryForms
    {             
      namespace Diffusion
      { 
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = 0.5 * int_u_v<Real, Scalar>(n, wt, u, v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = 0.5 * int_y_u_v<Real, Scalar>(n, wt, u, v, e);
          else 
            result = 0.5 * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          
          return result;
        }
        template
        scalar VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Jacobian::matrix_form(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;                                                               
        
        template<typename Real, typename Scalar>
        Scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          if (geom_type == HERMES_PLANAR) 
            result = 0.5 * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
          else if (geom_type == HERMES_AXISYM_X) 
            result = 0.5 * int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
          else 
            result = 0.5 * int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v, e);
          
          return result;
        }
        template
        scalar VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<scalar> *u_ext[],
                                                              Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;
        template
        Ord VacuumBoundaryCondition::Residual::vector_form(int n, double *wt, Func<Ord> *u_ext[],
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;        
       
        template<typename Real, typename Scalar>
        Scalar DiffusionReaction::Jacobian::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        {
          Scalar result;
          
          std::string mat = get_material(e->elem_marker, wf);     
          rank1 D_elem = matprop.get_D(mat);
          rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
          
          if (geom_type == HERMES_PLANAR) 
          {
            result = D_elem[g] * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) +
                      Sigma_r_elem[g] * int_u_v<Real, Scalar>(n, wt, u, v);
          }
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
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
        
        template<typename Real, typename Scalar>
        Scalar DiffusionReaction::Residual::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          Scalar result;
          
          std::string mat = get_material(e->elem_marker, wf);        
          rank1 D_elem = matprop.get_D(mat);
          rank1 Sigma_r_elem = matprop.get_Sigma_r(mat);
          
          if (geom_type == HERMES_PLANAR) 
          {
            result = D_elem[g] * int_grad_u_ext_grad_v<Real, Scalar>(n, wt, u_ext[g], v) +
                      Sigma_r_elem[g] * int_u_ext_v<Real, Scalar>(n, wt, u_ext[g], v);
          }
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
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
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Jacobian::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
        {
          if (!matprop.get_fission_multigroup_structure()[gto])
            return 0.0;
          
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          }
          
          std::string mat = get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::OuterIterationForm::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          if (!matprop.get_fission_multigroup_structure()[g])
            return 0.0;
            
          std::string mat = get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          if ((unsigned)ext->nf != nu_elem.size() || (unsigned)ext->nf != Sigma_f_elem.size())
            error(E_INVALID_GROUP_INDEX);
          
          Scalar result = 0;
          for (int i = 0; i < n; i++) 
          {
            Scalar local_res = 0;
            for (int gfrom = 0; gfrom < ext->nf; gfrom++)
              local_res += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext->fn[gfrom]->val[i];
            
            local_res = local_res * wt[i] * v->val[i];
            
            if (geom_type == HERMES_AXISYM_X)
              local_res = local_res * e->y[i];
            else if (geom_type == HERMES_AXISYM_Y)
              local_res = local_res * e->x[i];
            
            result += local_res;
          }
        
          return result * chi_elem[g] / keff;
        }
        
        template<typename Real, typename Scalar>
        Scalar FissionYield::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          if (!matprop.get_fission_multigroup_structure()[gto])
            return 0.0;
          
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
            else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
          }
          
          std::string mat = get_material(e->elem_marker, wf);
          rank1 nu_elem = matprop.get_nu(mat);
          rank1 Sigma_f_elem = matprop.get_Sigma_f(mat);
          rank1 chi_elem = matprop.get_chi(mat);
          
          return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar Scattering::Jacobian::matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const  
        {
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
            else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
          }
          
          return result * matprop.get_Sigma_s(get_material(e->elem_marker, wf))[gto][gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar Scattering::Residual::vector_form( int n, double *wt, Func<Scalar> *u_ext[],
                                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext ) const 
        { 
          Scalar result = 0;
          if (geom_type == HERMES_PLANAR) result = int_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) result = int_y_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
            else result = int_x_u_ext_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
          }
          
          return result * matprop.get_Sigma_s(get_material(e->elem_marker, wf))[gto][gfrom];
        }
        
        template<typename Real, typename Scalar>
        Scalar ExternalSources::LinearForm::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
        { 
          std::string mat = get_material(e->elem_marker, wf);
          
          if (geom_type == HERMES_PLANAR) 
            return matprop.get_src(mat)[g] * int_v<Real>(n, wt, v);
          else 
          {
            if (geom_type == HERMES_AXISYM_X) 
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
        void DefaultWeakFormFixedSource::lhs_init(unsigned int G, const MaterialPropertyMaps& matprop, 
                                                  GeomType geom_type)
        {
          bool2 Ss_nnz = matprop.get_scattering_multigroup_structure();
          bool1 chi_nnz = matprop.get_fission_multigroup_structure();
          
          for (unsigned int gto = 0; gto < G; gto++)
          {
            add_matrix_form(new DiffusionReaction::Jacobian(gto, matprop, geom_type));
            add_vector_form(new DiffusionReaction::Residual(gto, matprop, geom_type));
            
            for (unsigned int gfrom = 0; gfrom < G; gfrom++)
            {
              if (Ss_nnz[gto][gfrom])
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, matprop, geom_type));
              }
              
              if (chi_nnz[gto])
              {
                add_matrix_form(new FissionYield::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new FissionYield::Residual(gto, gfrom, matprop, geom_type));
              }
            }
          }
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                                               GeomType geom_type) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new ExternalSources::LinearForm(gto, matprop, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                DefaultFunction *f_src, std::string src_area,
                                                                GeomType geom_type  ) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, -1.0, f_src, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                DefaultFunction *f_src,
                                                                Hermes::vector<std::string> src_areas,
                                                                GeomType geom_type  ) : WeakForm(matprop.get_G())
        {
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, -1.0, f_src, geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                const std::vector<DefaultFunction*>& f_src,
                                                                std::string src_area, 
                                                                GeomType geom_type ) : WeakForm(matprop.get_G())
        {
          if (f_src.size() != matprop.get_G())
            error(E_INVALID_SIZE);
          
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_area, -1.0, f_src[gto], geom_type));
        }
        
        DefaultWeakFormFixedSource::DefaultWeakFormFixedSource( const MaterialPropertyMaps& matprop, 
                                                                const std::vector<DefaultFunction*>& f_src,
                                                                Hermes::vector<std::string> src_areas,
                                                                GeomType geom_type ) : WeakForm(matprop.get_G())
        {
          if (f_src.size() != matprop.get_G())
            error(E_INVALID_SIZE);
          
          lhs_init(matprop.get_G(), matprop, geom_type);
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            add_vector_form(new WeakFormsH1::DefaultVectorFormVol(gto, src_areas, -1.0, f_src[gto], geom_type));
        }
        
        DefaultWeakFormSourceIteration::DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                                                        Hermes::vector<MeshFunction*>& iterates,
                                                                        double initial_keff_guess, 
                                                                        GeomType geom_type ) : WeakForm(matprop.get_G())
        {      
          bool2 Ss_nnz = matprop.get_scattering_multigroup_structure();
          
          for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
          {
            add_matrix_form(new DiffusionReaction::Jacobian(gto, matprop, geom_type));
            add_vector_form(new DiffusionReaction::Residual(gto, matprop, geom_type));
            
            for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
            {
              if (Ss_nnz[gto][gfrom])
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, matprop, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, matprop, geom_type));
              }
            }
            
            FissionYield::OuterIterationForm* keff_iteration_form = 
              new FissionYield::OuterIterationForm( gto, matprop, iterates, initial_keff_guess, geom_type );
            keff_iteration_forms.push_back(keff_iteration_form);
            add_vector_form(keff_iteration_form);
          }
        }
        
        void DefaultWeakFormSourceIteration::update_keff(double new_keff) 
        { 
          std::vector<FissionYield::OuterIterationForm*>::iterator it = keff_iteration_forms.begin();
          for ( ; it != keff_iteration_forms.end(); ++it)
            (*it)->update_keff(new_keff); 
        }
      }
    }
  }
}
