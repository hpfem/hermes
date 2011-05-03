#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "h1.h"
#include "../function/forms.h"
#include <algorithm>

namespace WeakFormsNeutronics
{
  namespace Monoenergetic
  {    
    namespace Diffusion 
    {
      /* 
        Simple monoenergetic neutron diffusion, with the following weak formulation within each
        homogeneous region:
      
            \int_{region} D \nabla\phi \cdot \nabla\psi d\bfx + \int_{region} \Sigma_a \phi\psi d\bfx
              = \int_{region} Q_{ext}\psi d\bfx
        
        where 
        
            D         ... diffusion coefficient, 
            \Sigma_a  ... absorption cross-section, 
            Q_{ext}   ... external neutron sources 
          
        are region-wise constant physical parameters of the problem. Each region has one entry in vector
        'regions', which is the marker used for all elements it is composed of (usually specified in the
        mesh file). A corresponding entry in the *_map arguments is the value of the particular physical 
        parameter for that marker.
        
        Dirichlet and/or zero Neumann BC are assumed - nonzero Neumann or Newton boundary conditions can 
        be enabled by creating a descendant and adding surface forms to it.
      */
      class DefaultWeakFormFixedSource : public WeakForm
      {        
        public:
          DefaultWeakFormFixedSource( Hermes::vector<std::string> regions, 
                                      Hermes::vector<double> D_map, 
                                      Hermes::vector<double> Sigma_a_map, 
                                      Hermes::vector<double> Q_map ) : WeakForm(1) 
          {
            using namespace WeakFormsH1;
            
            for (unsigned int i = 0; i < regions.size(); i++)
            {
              /* Jacobian */
              // Diffusion.
              add_matrix_form(new DefaultJacobianDiffusion(0, 0, regions[i], D_map[i], HERMES_DEFAULT_SPLINE, HERMES_SYM));
              // Absorption.
              add_matrix_form(new DefaultMatrixFormVol(0, 0, regions[i], Sigma_a_map[i], HERMES_DEFAULT_FUNCTION, HERMES_SYM));
              
              /* Residual */
              // Diffusion.
              add_vector_form(new DefaultResidualDiffusion(0, regions[i], D_map[i]));
              // Absorption.
              add_vector_form(new DefaultResidualVol(0, regions[i], Sigma_a_map[i]));
              // Sources.
              add_vector_form(new DefaultVectorFormVol(0, regions[i], -Q_map[i]));
            }
          }
      };
    }
  }
    
  namespace Multigroup
  { 
    namespace MaterialProperties
    {
      namespace Definitions
      {
        typedef double rank0;
        typedef std::vector<double> rank1;
        typedef std::vector<std::vector<double > > rank2;
        typedef std::vector<std::vector<std::vector<double > > > rank3;
        
        typedef std::map<std::string, rank0> MaterialPropertyMap0;
        typedef std::map<std::string, rank1> MaterialPropertyMap1;
        typedef std::map<std::string, rank2> MaterialPropertyMap2;
        typedef std::map<std::string, rank3> MaterialPropertyMap3;
        
        typedef std::map<int,rank0> iMarkerPropertyMap0;
        typedef std::map<int,rank1> iMarkerPropertyMap1;
        typedef std::map<int,rank2> iMarkerPropertyMap2;
        typedef std::map<int,rank3> iMarkerPropertyMap3;
      }
      
      namespace Messages
      {
        static const char* E_INF_VALUE = 
          "Attempt to set an infinite material property.";
        static const char* W_NEG_VALUE =
          "Entered material data lead to some negative properties.";
        static const char* W_MG_EXTENSION = 
          "Attempted to create a multigroup material-property map in a container for singlegroup maps.";
        static const char* W_SA_LT_SF =
          "Possible unphysical situation detected: Sigma_a < Sigma_f.";  
        static const char* E_MR_EXTENSION = 
          "Cannot create a multiregion material-property map: no regions specified.";
        static const char* E_INSUFFICIENT_DATA =
          "Not all required material properties have been set.";
        static const char* W_NO_FISSION =
          "Not all required fission properties have been set or could be determined automatically."
          "Assuming a non-fissioning system.";
        static const char* W_NO_SCATTERING =
          "Not all required scattering properties have been set or could be determined automatically."
          "Assuming a purely absorbing system.";
        static const char* E_INVALID_COMBINATION =
          "Invalid combination of entered material properties.";
        static const char* E_NONMATCHING_PROPERTIES =
          "All properties must be defined for a single given number of materials.";
        static const char* E_INVALID_SIZE =
          "Material property defined for an unexpected number of groups.";
        static const char* E_INVALID_GROUP_INDEX =
          "Attempted to access an out-of-range group.";
        static const char* E_INVALID_MARKER =
          "Material data undefined for the given element marker.";
        static const char* E_SG_SIGMA_R = 
          "Group-reduction cross-section (Sigma_r) is not defined for one-group (i.e. monoenergetic) problems."
          "Set Sigma_a instead.";
      }
      
      namespace ValidationFunctors
      {
        using namespace Definitions;
        using namespace Messages;
        
        struct ensure_trivial { 
          void operator() (MaterialPropertyMap1::value_type x) { 
            MaterialPropertyMap1::mapped_type::iterator it;
            for (it = x.second.begin(); it != x.second.end(); ++it) 
              if (fabs(*it) > 1e-14)
                error(E_INVALID_COMBINATION);
          }
        };
        
        struct ensure_size { 
          ensure_size(unsigned int nrows, unsigned int ncols = 0) 
            : nrows(nrows), ncols(ncols) {};
          
          void operator() (MaterialPropertyMap1::value_type x) { 
            if (x.second.size() != nrows)
              error(E_INVALID_SIZE);
          }
          
          void operator() (MaterialPropertyMap2::value_type x) {
            if (x.second.size() != nrows)
              error(E_INVALID_SIZE);
            
            MaterialPropertyMap2::mapped_type::iterator it;
            for (it = x.second.begin(); it != x.second.end(); ++it) 
              if (it->size() != ncols)
                error(E_INVALID_SIZE);
          }
          
          private:
            unsigned int nrows, ncols;
        };
        
        struct ensure_size_at_least { 
          ensure_size_at_least(unsigned int nrows, unsigned int ncols = 0) 
            : nrows(nrows), ncols(ncols) {};
          
          void operator() (iMarkerPropertyMap1::value_type x) { 
            if (x.second.size() < nrows)
              error(E_INVALID_GROUP_INDEX);
          }
          
          void operator() (iMarkerPropertyMap2::value_type x) { 
            if (x.second.size() < nrows || x.second[nrows].size() < ncols)
              error(E_INVALID_GROUP_INDEX);
          }
          
          private:
            unsigned int nrows, ncols;
        };
      }
      
      namespace Common
      {
        using namespace Definitions;
        using namespace Messages;
        
        class NDArrayMapOp
        {
          //
          // NOTE: Could be perhaps combined with the classes material_property_map and MultiArray below
          // and moved to hermes_common as a general way of handling maps with multidimensional mapped types.
          //
          
          template <typename NDArrayType>
          static rank0 divide(rank0 x, rank0 y) 
          {
            if (x == 0 && y == 0) 
              return 0.0;
            else if (y == 0)
            {
              error(E_INF_VALUE);
              return -1.0;
            }
            else
              return x/y;
          }
          
          template <typename NDArrayType>
          static rank0 multiply(rank0 x, rank0 y) 
          {
            return x*y;
          }
          
          template <typename NDArrayType>
          static rank0 add(rank0 x, rank0 y) 
          {
            return x + y;
          }
          
          template <typename NDArrayType>
          static rank0 subtract(rank0 x, rank0 y) 
          {
            rank0 ret = x - y;
            if(ret < 0)
              warning(W_NEG_VALUE);
            return ret;
          }
          
          #define for_each_element_in_dimension \
                      typedef typename NDArrayType::value_type dim_type;                      \
                      typename NDArrayType::const_iterator dim_iterator_x = x.begin();        \
                      typename NDArrayType::const_iterator dim_iterator_y = y.begin();        \
                      for ( ; dim_iterator_x != x.end(); ++dim_iterator_x, ++dim_iterator_y )                    
          
          template <typename NDArrayType>
          static NDArrayType divide(const NDArrayType& x, const NDArrayType& y) 
          { 
            NDArrayType res; res.reserve(x.size());
                  
            for_each_element_in_dimension
              res.push_back( divide<dim_type>(*dim_iterator_x, *dim_iterator_y) );
            
            return res;
          }
          
          template <typename NDArrayType>
          static NDArrayType multiply(const NDArrayType& x, const NDArrayType& y) 
          { 
            NDArrayType res; res.reserve(x.size());
            
            for_each_element_in_dimension
              res.push_back( multiply<dim_type>(*dim_iterator_x, *dim_iterator_y) );
            
            return res;
          }
          
          template <typename NDArrayType>
          static NDArrayType add(const NDArrayType& x, const NDArrayType& y) 
          { 
            NDArrayType res; res.reserve(x.size());
            
            for_each_element_in_dimension
              res.push_back( add<dim_type>(*dim_iterator_x, *dim_iterator_y) );
            
            return res;
          }
          
          template <typename NDArrayType>
          static NDArrayType subtract(const NDArrayType& x, const NDArrayType& y) 
          { 
            NDArrayType res; res.reserve(x.size());
            
            for_each_element_in_dimension
              res.push_back( subtract<dim_type>(*dim_iterator_x, *dim_iterator_y) );
            
            return res;
          }
          
          
          #undef for_each_element_in_dimension
          
          #define for_each_element_in_map \
                    typename std::map<std::string, T>::iterator iterator_ret = ret.begin();   \
                    typename std::map<std::string, T>::const_iterator iterator_x = x.begin(); \
                    typename std::map<std::string, T>::const_iterator iterator_y = y.begin(); \
                    for ( ; iterator_x != x.end(); ++iterator_x, ++iterator_y, ++iterator_ret )                    
          
          public:
            template <typename T>
            static std::map<std::string, T> divide(const std::map<std::string, T>& x, 
                                                   const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = divide<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
            
            template <typename T>
            static std::map<std::string, T> multiply(const std::map<std::string, T>& x, 
                                                     const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = multiply<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
                        
            template <typename T>
            static std::map<std::string, T> add(const std::map<std::string, T>& x, 
                                                const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = add<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
            
            template <typename T>
            static std::map<std::string, T> subtract(const std::map<std::string, T>& x, 
                                                     const std::map<std::string, T>& y)
            {
              std::map<std::string, T> ret = x;
              
              for_each_element_in_map
                iterator_ret->second = subtract<T>(iterator_x->second, iterator_y->second);
              
              return ret;
            }
            
            #undef for_each_element_in_map
        };
        
        
        class MaterialPropertyMaps
        {
          protected:
                                
            MaterialPropertyMap1 Sigma_f;
            MaterialPropertyMap1 nu;
            MaterialPropertyMap1 chi;
            
            MaterialPropertyMap1 Sigma_a;
            MaterialPropertyMap1 nuSigma_f;
            
            std::set<std::string> materials_list;
            unsigned int G;
                  
            void extend_to_multigroup(const MaterialPropertyMap0& mrsg_map, MaterialPropertyMap1 *mrmg_map)
            {
              if (G == 1)
                warning(W_MG_EXTENSION);
              
              MaterialPropertyMap0::const_iterator it;
              for (it = mrsg_map.begin(); it != mrsg_map.end(); ++it)
                mrmg_map->at(it->first).assign(G, it->second);
            
            }
            
            void extend_to_multiregion(const rank1& srmg_array, MaterialPropertyMap1 *mrmg_map)
            {
              if (materials_list.empty())
                error(E_MR_EXTENSION);
              
              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                mrmg_map->at(*it) = srmg_array;
            }
            
            void extend_to_multiregion_multigroup(const rank0& srsg_value, MaterialPropertyMap1 *mrmg_map)
            {
              if (materials_list.empty())
                error(E_MR_EXTENSION);
              if (G == 1)
                warning(W_MG_EXTENSION);
              
              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                mrmg_map->at(*it).assign(G, srsg_value);
            }
            
            void fill_with(double c, MaterialPropertyMap1 *mrmg_map)
            {
              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                mrmg_map->at(*it).assign(G, c);
            }
                      
            MaterialPropertyMaps(unsigned int G, std::set<std::string> mat_list = std::set<std::string>()) 
              : materials_list(mat_list), G(G)  { };
                        
            virtual void validate()
            {       
              using namespace ValidationFunctors;
              
              if (chi.empty())
                fill_with(1.0, &chi);
              
              if (nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
                nu = NDArrayMapOp::divide<rank1>(nuSigma_f, Sigma_f);
              else if (nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
                nuSigma_f = NDArrayMapOp::multiply<rank1>(nu, Sigma_f);
              else if (Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
                Sigma_f = NDArrayMapOp::divide<rank1>(nuSigma_f, nu);
              else if (!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
              {
                MaterialPropertyMap1 diff = NDArrayMapOp::subtract<rank1>( nuSigma_f, 
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
            
          public:
            
            void set_nu(const MaterialPropertyMap1& nu)
            {
              this->nu = nu;
            }
            
            void set_nu(const MaterialPropertyMap0& nu)
            {
              extend_to_multigroup(nu, &this->nu);      
            }
            
            void set_nu(const rank1& nu)
            {
              extend_to_multiregion(nu, &this->nu);
            }
            
            void set_nu(const rank0& nu)
            {
              extend_to_multiregion_multigroup(nu, &this->nu);
            }
                  
            void set_chi(const MaterialPropertyMap1& chi)
            {
              this->chi = chi;
            }
            
            void set_chi(const rank1& chi)
            {
              extend_to_multiregion(chi, &this->chi);
            }
            
            void set_Sigma_a(const MaterialPropertyMap1& Sa)
            {
              this->Sigma_a = Sa;
            }
            
            void set_Sigma_f(const MaterialPropertyMap1& Sf)
            {
              this->Sigma_f = Sf;
            }
            
            void set_nuSigma_f(const MaterialPropertyMap1 nSf)
            {
              this->nuSigma_f = nSf;
            }
            
            const MaterialPropertyMap1& get_Sigma_f() const
            {
              return this->Sigma_f;
            }
            const MaterialPropertyMap1& get_nu() const
            {
              return this->nu;
            }
            const MaterialPropertyMap1& get_chi() const
            {
              return this->chi;
            }
            
            unsigned int get_G() const { return G; }          
        };
      }
      
      namespace Diffusion
      {
        using namespace Definitions;
        using namespace Messages;
        
        class MaterialPropertyMaps : public Common::MaterialPropertyMaps
        {
          protected:
            
            MaterialPropertyMap1 D;
            MaterialPropertyMap1 Sigma_r;
            MaterialPropertyMap2 Sigma_s;
            MaterialPropertyMap1 src;
            
            MaterialPropertyMap1 Sigma_t;
            
          public:
            
            MaterialPropertyMaps(unsigned int G, std::set<std::string> mat_list = std::set<std::string>()) 
              : Common::MaterialPropertyMaps(G, mat_list) { };
              
            MaterialPropertyMap1 extract_map2_diagonals(const MaterialPropertyMap2& map2)
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
            
            MaterialPropertyMap1 sum_map2_columns(const MaterialPropertyMap2& map2)
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
            
            MaterialPropertyMap2 create_map2_by_diagonals(const MaterialPropertyMap1& diags)
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
            
            void fill_with(double c, MaterialPropertyMap2 *mrmg_map)
            {
              std::set<std::string>::const_iterator it;
              for (it = materials_list.begin(); it != materials_list.end(); ++it)
                mrmg_map->at(*it).assign(G, rank1(G, c));
            }
            
            // We always need to supply chi, nu, Sigma_f, Sigma_r, Sigma_s and D to our neutronics weak forms. 
            // These parameters are often defined in terms of the other ones, or not specified at all and assumed 
            // to be zero for a particular simplified situation. This method, together with its complement in the
            // parent class, uses the most typical definitions to build the six-parameter set from the given input. 
            // It also checks whether the user did not enter nonsensical values. However, values entered by the 
            // user may sometimes not satisfy the common relations, as some empirical corrections may have been 
            // already included in them.
            virtual void validate()
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
              
              if (!Sigma_s_given)
              {
                // If Sigma_s is not given, but Sigma_t is, we can obtain the former from the latter and from Sigma_r.
                // Note that the execution will come here only if the user entered Sigma_r himself - otherwise, Sigma_s
                // has been already set in the previous test case.
                
                if (Sigma_t_given)
                  Sigma_s = create_map2_by_diagonals(Common::NDArrayMapOp::subtract<rank1>(Sigma_t, Sigma_r));
                else
                {
                  warning(W_NO_SCATTERING);
                  fill_with(0.0, &Sigma_s);
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
            
            void set_src(const MaterialPropertyMap1& src)
            {
              this->src = src;
            }
            
            void set_src(const MaterialPropertyMap0& src)
            {
              extend_to_multigroup(src, &this->src);            
            }
            
            void set_src(const rank1& src)
            {
              extend_to_multiregion(src, &this->src);
            }
            
            void set_src(const double& src)
            {
              extend_to_multiregion_multigroup(src, &this->src);
            }
            
            void set_D(const MaterialPropertyMap1& D)
            {
              this->D = D;
            }
            
            void set_Sigma_r(const MaterialPropertyMap1& Sr)
            {
              this->Sigma_r = Sr;
            }
            
            void set_Sigma_t(const MaterialPropertyMap1& St)
            {
              this->Sigma_t = St;
            }
            
            void set_Sigma_s(const MaterialPropertyMap2& Ss)
            {
              this->Sigma_s = Ss;
            }
            
            const MaterialPropertyMap2& get_Sigma_s() const
            {
              return this->Sigma_s;
            }
            const MaterialPropertyMap1& get_Sigma_r() const
            {
              return this->Sigma_r;
            }
            const MaterialPropertyMap1& get_D() const
            {
              return this->D;
            }
            const MaterialPropertyMap1& get_src() const
            {
              return this->src;
            }
        };
      }  
    
      template <typename NDArrayType>
      class material_property_map
      {
        private:
          std::map<std::string, NDArrayType> m_map;
        public:
          material_property_map(const std::string& key, const NDArrayType& val)
          {
            m_map[key] = val;
          }
          
          material_property_map<NDArrayType>& operator()(const std::string& key, const NDArrayType& val)
          {
            m_map[key] = val;
            return *this;
          }
          
          operator std::map<std::string, NDArrayType>()
          {
            return m_map;
          }
      };
      
      template <typename NDArrayType>
      class MultiArray
      {
        private:
          std::vector<NDArrayType> m_data;
        public:
          MultiArray(const NDArrayType& val)
          {
            m_data.push_back(val);
          }
          
          MultiArray<NDArrayType>& operator()(const NDArrayType& val)
          {
            m_data.push_back(val);
            return *this;
          }
          
          operator std::vector<NDArrayType>()
          {
            return m_data;
          }
      };
      
      namespace Definitions
      {
        typedef MultiArray<rank0> grow;
        typedef MultiArray<rank1> gmat;
      }
    }
                                 
    namespace ElementaryForms
    {
      using namespace MaterialProperties::Definitions;
      
      namespace Diffusion
      { 
        using namespace MaterialProperties::Diffusion;
        
        class DiffusionReaction
        {
          public:
            
            class Data
            {
              iMarkerPropertyMap1 D;
              iMarkerPropertyMap1 Sigma_r;
              
              public:
                Data() {};
                
                Data(const MaterialPropertyMaps& matprop, 
                     Mesh::ElementMarkersConversion *element_markers_conversion) 
                {
                  MaterialPropertyMap1 str_mpm = matprop.get_D();
                  MaterialPropertyMap1::const_iterator str_mpm_iter = str_mpm.begin();
                  for ( ; str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    D[int_marker] = str_mpm_iter->second;
                  }
                  
                  str_mpm = matprop.get_Sigma_r();
                  for (str_mpm_iter = str_mpm.begin(); str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    Sigma_r[int_marker] = str_mpm_iter->second;
                  }
                  
                  if (D.size() != Sigma_r.size()) 
                    error(E_NONMATCHING_PROPERTIES);
                }
                
                void check_validity(unsigned int g) const
                {
                  using MaterialProperties::ValidationFunctors::ensure_size_at_least;
                  
                  std::for_each(D.begin(), D.end(), ensure_size_at_least(g));
                  std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size_at_least(g));
                }
                
                const rank1& D_elem(int elem_marker) const
                {
                  // Constant properties within the current active element. Note that prop[e->elem_marker]
                  // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
                  iMarkerPropertyMap1::const_iterator data = D.find(elem_marker);
                  if (data != D.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
                
                const rank1& Sigma_r_elem(int elem_marker) const
                {
                  iMarkerPropertyMap1::const_iterator data = Sigma_r.find(elem_marker);
                  if (data != Sigma_r.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
            };
            
            class Jacobian : public WeakForm::MatrixFormVol
            {
              public:            
                Jacobian(unsigned int g, 
                        const Data& data, GeomType geom_type = HERMES_PLANAR) 
                  : WeakForm::MatrixFormVol(g, g, HERMES_ANY, HERMES_SYM),
                    g(g), data(data), geom_type(geom_type)
                {
                  data.check_validity(g);
                }
                Jacobian(unsigned int g, std::string area,
                        const Data& data, GeomType geom_type = HERMES_PLANAR)
                  : WeakForm::MatrixFormVol(g, g, area, HERMES_SYM),
                    g(g), data(data), geom_type(geom_type)
                {  
                  data.check_validity(g);
                }
                
                template<typename Real, typename Scalar>
                Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
                {
                  Scalar result = 0;
                  
                  rank1 D_elem = data.D_elem(e->elem_marker);
                  rank1 Sigma_r_elem = data.Sigma_r_elem(e->elem_marker);
                  
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

                virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                      Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const 
                { 
                  return  matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
                }
                
                virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const 
                { 
                  return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
                }

                // This is to make the form usable in rk_time_step().
                virtual WeakForm::MatrixFormVol* clone() {
                  return new Jacobian(*this);
                }

              private:
                
                unsigned int g;
                const Data& data;
                GeomType geom_type;
            };
            
            class Residual : public WeakForm::VectorFormVol
            {
              public:
                
                Residual(unsigned int g, 
                         const Data& data, GeomType geom_type = HERMES_PLANAR) 
                  : WeakForm::VectorFormVol(g, HERMES_ANY),
                    g(g), data(data), geom_type(geom_type)
                {
                  data.check_validity(g);
                }
                Residual(unsigned int g, std::string area,
                         const Data& data, GeomType geom_type = HERMES_PLANAR)
                  : WeakForm::VectorFormVol(g, area),
                    g(g), data(data), geom_type(geom_type)
                {  
                  data.check_validity(g);
                }
                
                template<typename Real, typename Scalar>
                Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
                { 
                  Scalar result = 0;
                  
                  rank1 D_elem = data.D_elem(e->elem_marker);
                  rank1 Sigma_r_elem = data.Sigma_r_elem(e->elem_marker);
                  
                  if (geom_type == HERMES_PLANAR) 
                  {
                    result = D_elem[g] * int_grad_u_grad_v<Real, Scalar>(n, wt, u_ext[g], v) +
                             Sigma_r_elem[g] * int_u_v<Real, Scalar>(n, wt, u_ext[g], v);
                  }
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) 
                    {
                      result = D_elem[g] * int_y_grad_u_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) + 
                               Sigma_r_elem[g] * int_y_u_v<Real, Scalar>(n, wt, u_ext[g], v, e);
                    }
                    else 
                    {
                      result = D_elem[g] * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u_ext[g], v, e) + 
                               Sigma_r_elem[g] * int_x_u_v<Real, Scalar>(n, wt, u_ext[g], v, e);
                    }
                  }
                  return result;
                }
                
                virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
                {
                  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
                {
                  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
                }
                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::VectorFormVol* clone() {
                  return new Residual(*this);
                }
                
              private:
                
                unsigned int g;
                const Data& data;
                GeomType geom_type;
            };
        };
      
        class FissionYield
        {
          public:
            
            class Data
            {
              iMarkerPropertyMap1 chi;
              iMarkerPropertyMap1 Sigma_f;
              iMarkerPropertyMap1 nu;
              
              public:
                Data() {};
                
                Data(const MaterialPropertyMaps& matprop,
                     Mesh::ElementMarkersConversion *element_markers_conversion)
                {
                  MaterialPropertyMap1 str_mpm = matprop.get_Sigma_f();
                  MaterialPropertyMap1::const_iterator str_mpm_iter = str_mpm.begin();
                  for ( ; str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    Sigma_f[int_marker] = str_mpm_iter->second;
                  }
                  
                  str_mpm = matprop.get_nu();
                  for (str_mpm_iter = str_mpm.begin(); str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    nu[int_marker] = str_mpm_iter->second;
                  }
                  
                  str_mpm = matprop.get_chi();
                  for (str_mpm_iter = str_mpm.begin(); str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    chi[int_marker] = str_mpm_iter->second;
                  }
                  
                  if (chi.size() != Sigma_f.size() || chi.size() != nu.size())
                    error(E_NONMATCHING_PROPERTIES);
                }
                
                void check_validity(unsigned int gto, unsigned int gfrom = 0) const
                {
                  using MaterialProperties::ValidationFunctors::ensure_size_at_least;
                  
                  std::for_each(chi.begin(), chi.end(), ensure_size_at_least(gto));
                  
                  if (gfrom > 0) {
                    std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size_at_least(gfrom));
                    std::for_each(nu.begin(), nu.end(), ensure_size_at_least(gfrom));
                  }
                }
                
                const rank1& chi_elem(int elem_marker) const
                {
                  // Constant properties within the current active element. Note that prop[e->elem_marker]
                  // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
                  iMarkerPropertyMap1::const_iterator data = chi.find(elem_marker);
                  if (data != chi.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
                
                const rank1& Sigma_f_elem(int elem_marker) const
                {
                  iMarkerPropertyMap1::const_iterator data = Sigma_f.find(elem_marker);
                  if (data != Sigma_f.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
                
                const rank1& nu_elem(int elem_marker) const
                {
                  iMarkerPropertyMap1::const_iterator data = nu.find(elem_marker);
                  if (data != nu.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
            };
          
            class Jacobian : public WeakForm::MatrixFormVol
            {
              public:
                
                Jacobian( unsigned int gto, unsigned int gfrom, 
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::MatrixFormVol(gto, gfrom), 
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                {
                  data.check_validity(gto, gfrom);
                }
                
                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::MatrixFormVol(gto, gfrom, area),
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                { 
                  data.check_validity(gto, gfrom);
                }
                
                template<typename Real, typename Scalar>
                Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
                {
                  Scalar result = 0;
                  if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
                    else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
                  }
                  
                  rank1 nu_elem = data.nu_elem(e->elem_marker);
                  rank1 Sigma_f_elem = data.Sigma_f_elem(e->elem_marker);
                  rank1 chi_elem = data.chi_elem(e->elem_marker);
                  
                  return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
                }
                
                virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                      Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const 
                { 
                  return  -1 * matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
                }
                
                virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const 
                { 
                  return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
                }
                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::MatrixFormVol* clone() {
                  return new Jacobian(*this);
                }
                
              private:
                
                unsigned int gto, gfrom;
                const Data& data;
                GeomType geom_type;
            };
        
            class OuterIterationForm : public WeakForm::VectorFormVol
            {
              public:
                
                OuterIterationForm( unsigned int g, 
                                    const Data& data, double keff = 1.0,
                                    GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(g),
                    g(g), keff(keff), data(data), geom_type(geom_type)
                {
                  data.check_validity(g);
                }
                
                OuterIterationForm( unsigned int g, std::string area,
                                    const Data& data, double keff = 1.0,
                                    GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(g, area),
                    g(g), keff(keff), data(data), geom_type(geom_type)
                {
                  data.check_validity(g);
                }
                
                template<typename Real, typename Scalar>
                Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
                { 
                  rank1 nu_elem = data.nu_elem(e->elem_marker);
                  rank1 Sigma_f_elem = data.Sigma_f_elem(e->elem_marker);
                  rank1 chi_elem = data.chi_elem(e->elem_marker);
                  
                  if ((unsigned)ext->nf != nu_elem.size() || (unsigned)ext->nf != Sigma_f_elem.size())
                    error(E_INVALID_GROUP_INDEX);
                  
                  Scalar result = 0;
                  for (int i = 0; i < n; i++) 
                  {
                    for (int gfrom = 0; gfrom < ext->nf; gfrom++)
                      result += nu_elem[gfrom] * Sigma_f_elem[gfrom] * ext->fn[gfrom]->val[i];
                    
                    result = result * wt[i] * v->val[i];
                    
                    if (geom_type == HERMES_AXISYM_X)
                      result = result * e->y[i];
                    else if (geom_type == HERMES_AXISYM_Y)
                      result = result * e->x[i];
                  }
                  
                  return result * chi_elem[g] / keff;
                }

                virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const 
                {
                  return -1 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const 
                {
                  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step().
                virtual WeakForm::VectorFormVol* clone() {
                  return new OuterIterationForm(*this);
                }
                
                void update_keff(double new_keff) { keff = new_keff; }
                
              private:
                
                unsigned int g;
                double keff;
                const Data& data;
                GeomType geom_type;
            };
          
            class Residual : public WeakForm::VectorFormVol
            {
              public:
                Residual( unsigned int gto, unsigned int gfrom, 
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(gto), 
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                {
                  data.check_validity(gto, gfrom);
                }
                
                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(gto, area),
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                { 
                  data.check_validity(gto, gfrom);
                }
                
                template<typename Real, typename Scalar>
                Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
                { 
                  Scalar result = 0;
                  if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                    else result = int_x_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                  }
                  
                  rank1 nu_elem = data.nu_elem(e->elem_marker);
                  rank1 Sigma_f_elem = data.Sigma_f_elem(e->elem_marker);
                  rank1 chi_elem = data.chi_elem(e->elem_marker);
                  
                  return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
                }
                
                virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
                {
                  return -1 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
                {
                  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
                }
                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::VectorFormVol* clone() {
                  return new Residual(*this);
                }
                
              private:
                
                unsigned int gto, gfrom;
                const Data& data;
                GeomType geom_type;
            };
        };
  
        class Scattering
        {
          public:
            
            class Data
            {
              iMarkerPropertyMap2 Sigma_s;
              
              public:
                Data() {};
                
                Data(const MaterialPropertyMaps& matprop,
                     Mesh::ElementMarkersConversion *element_markers_conversion) 
                {
                  MaterialPropertyMap2 str_mpm2 = matprop.get_Sigma_s();
                  MaterialPropertyMap2::const_iterator str_mpm2_iter = str_mpm2.begin();
                  for ( ; str_mpm2_iter != str_mpm2.end(); ++str_mpm2_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm2_iter->first);
                    Sigma_s[int_marker] = str_mpm2_iter->second;
                  }
                }
                
                void check_validity(unsigned int gto, unsigned int gfrom) const
                {
                  std::for_each(Sigma_s.begin(), Sigma_s.end(), 
                                MaterialProperties::ValidationFunctors::ensure_size_at_least(gto, gfrom));
                }
                
                const rank2& Sigma_s_elem(int elem_marker) const
                {
                  // Constant properties within the current active element. Note that prop[e->elem_marker]
                  // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
                  iMarkerPropertyMap2::const_iterator data = Sigma_s.find(elem_marker);
                  if (data != Sigma_s.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank2()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
            };
        
            class Jacobian : public WeakForm::MatrixFormVol
            {
              public:
                
                Jacobian( unsigned int gto, unsigned int gfrom, 
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::MatrixFormVol(gto, gfrom), 
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                {
                  data.check_validity(gto, gfrom);
                }
                
                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::MatrixFormVol(gto, gfrom, area),
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                { 
                  data.check_validity(gto, gfrom);
                }
                
                template<typename Real, typename Scalar>
                Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const  
                {
                  Scalar result = 0;
                  if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u, v);
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u, v, e);
                    else result = int_x_u_v<Real, Scalar>(n, wt, u, v, e);
                  }
                 
                  return result * data.Sigma_s_elem(e->elem_marker)[gto][gfrom];
                }
                
                virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                      Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const 
                { 
                  return  -1 * matrix_form<double, scalar> (n, wt, u_ext, u, v, e, ext);
                }
                
                virtual Ord ord( int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext  ) const 
                { 
                  return  matrix_form<Ord, Ord> (n, wt, u_ext, u, v, e, ext);
                }
                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::MatrixFormVol* clone() {
                  return new Jacobian(*this);
                }
                
              private:
                
                unsigned int gto, gfrom;
                const Data& data;
                GeomType geom_type;
            };
          
            class Residual : public WeakForm::VectorFormVol
            {
              public:
                Residual( unsigned int gto, unsigned int gfrom, 
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(gto), 
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                {
                  data.check_validity(gto, gfrom);
                }
                
                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                          const Data& data, GeomType geom_type = HERMES_PLANAR )
                  : WeakForm::VectorFormVol(gto, area),
                    gto(gto), gfrom(gfrom), data(data), geom_type(geom_type)
                { 
                  data.check_validity(gto, gfrom);
                }
                
                template<typename Real, typename Scalar>
                Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
                { 
                  Scalar result = 0;
                  if (geom_type == HERMES_PLANAR) result = int_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v);
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                    else result = int_x_u_v<Real, Scalar>(n, wt, u_ext[gfrom], v, e);
                  }
                  
                  return result * data.Sigma_s_elem(e->elem_marker)[gto][gfrom];
                }
                
                virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
                {
                  return -1 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
                {
                  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
                }
                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::VectorFormVol* clone() {
                  return new Residual(*this);
                }
                
              private:
                
                unsigned int gto, gfrom;
                const Data& data;
                GeomType geom_type;
            };
        };
        
        class ExternalSources
        {
          public:
            
            class Data
            {
              iMarkerPropertyMap1 src;
              
              public:
                Data() {};
                
                Data(const MaterialPropertyMaps& matprop,
                     Mesh::ElementMarkersConversion *element_markers_conversion) 
                {
                  MaterialPropertyMap1 str_mpm = matprop.get_src();
                  MaterialPropertyMap1::const_iterator str_mpm_iter = str_mpm.begin();
                  for ( ; str_mpm_iter != str_mpm.end(); ++str_mpm_iter)
                  {
                    int int_marker = element_markers_conversion->get_internal_marker(str_mpm_iter->first);
                    src[int_marker] = str_mpm_iter->second;
                  }
                }
                
                void check_validity(unsigned int g) const
                {
                  std::for_each(src.begin(), src.end(), 
                                MaterialProperties::ValidationFunctors::ensure_size_at_least(g));
                }
                
                const rank1& src_elem(int elem_marker) const
                {
                  // Constant properties within the current active element. Note that prop[e->elem_marker]
                  // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
                  iMarkerPropertyMap1::const_iterator data = src.find(elem_marker);
                  if (data != src.end())
                    return data->second;
                  else
                  {
                    error(E_INVALID_MARKER);
                    return *(new rank1()); // To avoid MSVC problems; execution should never come to this point.
                  }
                }
            };
        
            class LinearForm : public WeakForm::VectorFormVol
            {
              public:
                
                LinearForm( unsigned int g, 
                            const Data& data, GeomType geom_type = HERMES_PLANAR)
                  : WeakForm::VectorFormVol(g), g(g), data(data), geom_type(geom_type) 
                { 
                  data.check_validity(g);
                }
                
                LinearForm( unsigned int g, std::string area,
                            const Data& data, GeomType geom_type = HERMES_PLANAR)
                  : WeakForm::VectorFormVol(g, area), g(g), data(data), geom_type(geom_type)
                { 
                  data.check_validity(g);
                }
                
                template<typename Real, typename Scalar>
                Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                  Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
                { 
                  if (geom_type == HERMES_PLANAR) 
                    return data.src_elem(e->elem_marker)[g] * int_v<Real>(n, wt, v);
                  else 
                  {
                    if (geom_type == HERMES_AXISYM_X) 
                      return data.src_elem(e->elem_marker)[g] * int_y_v<Real>(n, wt, v, e);
                    else 
                      return data.src_elem(e->elem_marker)[g] * int_x_v<Real>(n, wt, v, e);
                  }
                }
                
                virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
                {
                  return -1 * vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
                {
                  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
                }
                                
                // This is to make the form usable in rk_time_step().
                virtual WeakForm::VectorFormVol* clone() {
                  return new LinearForm(*this);
                }
              
              private:
                
                unsigned int g;
                const Data& data;
                GeomType geom_type;        
            }; 
        };
            
      }                
    
    }
    
    namespace CompleteWeakForms
    { 
      using namespace MaterialProperties::Definitions;
            
      namespace Diffusion
      {      
        using namespace MaterialProperties::Diffusion;
        using namespace ElementaryForms::Diffusion;
        
        class GenericMultigroupDiffusionWeakForm : public WeakForm
        {
          protected:
            
            DiffusionReaction::Data diffusion_reaction_data;
            FissionYield::Data fission_yield_data;
            Scattering::Data scattering_data;
            ExternalSources::Data ext_src_data;
            
          public:
            
            GenericMultigroupDiffusionWeakForm(const MaterialPropertyMaps& matprop) 
              : diffusion_reaction_data(matprop, element_markers_conversion),
                fission_yield_data(matprop, element_markers_conversion),
                scattering_data(matprop, element_markers_conversion),
                ext_src_data(matprop, element_markers_conversion)
            {};
        };
        
        class DefaultWeakFormFixedSource : public GenericMultigroupDiffusionWeakForm
        { 
          DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                     GeomType geom_type = HERMES_PLANAR) 
            : GenericMultigroupDiffusionWeakForm(matprop)
          {            
            for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            {
              add_matrix_form(new DiffusionReaction::Jacobian(gto, diffusion_reaction_data, geom_type));
              add_vector_form(new DiffusionReaction::Residual(gto, diffusion_reaction_data, geom_type));
              
              for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, scattering_data, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, scattering_data, geom_type));
                
                add_matrix_form(new FissionYield::Jacobian(gto, gfrom, fission_yield_data, geom_type));
                add_vector_form(new FissionYield::Residual(gto, gfrom, fission_yield_data, geom_type));
              }
              
              add_vector_form(new ExternalSources::LinearForm(gto, ext_src_data, geom_type));
            }
          }
          
          DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, 
                                     DefaultFunction *f_src,
                                     GeomType geom_type = HERMES_PLANAR) 
            : GenericMultigroupDiffusionWeakForm(matprop)
          {            
            for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            {
              add_matrix_form(new DiffusionReaction::Jacobian(gto, diffusion_reaction_data, geom_type));
              add_vector_form(new DiffusionReaction::Residual(gto, diffusion_reaction_data, geom_type));
              
              for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, scattering_data, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, scattering_data, geom_type));
                
                add_matrix_form(new FissionYield::Jacobian(gto, gfrom, fission_yield_data, geom_type));
                add_vector_form(new FissionYield::Residual(gto, gfrom, fission_yield_data, geom_type));
              }
              
              add_vector_form(new DefaultVectorFormVol(gto, 1.0, f_src, geom_type));
            }
          }
        };
                
        class DefaultWeakFormSourceIteration : public GenericMultigroupDiffusionWeakForm
        {
          DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop,
                                          double initial_keff_guess,
                                          GeomType geom_type = HERMES_PLANAR ) 
            : GenericMultigroupDiffusionWeakForm(matprop)
          {            
            for (unsigned int gto = 0; gto < matprop.get_G(); gto++)
            {
              add_matrix_form(new DiffusionReaction::Jacobian(gto, diffusion_reaction_data, geom_type));
              add_vector_form(new DiffusionReaction::Residual(gto, diffusion_reaction_data, geom_type));
              
              for (unsigned int gfrom = 0; gfrom < matprop.get_G(); gfrom++)
              {
                add_matrix_form(new Scattering::Jacobian(gto, gfrom, scattering_data, geom_type));
                add_vector_form(new Scattering::Residual(gto, gfrom, scattering_data, geom_type));
              }
              
              add_vector_form(new FissionYield::OuterIterationForm(gto, fission_yield_data, 
                                                                   initial_keff_guess, geom_type));
            }
          }
        };
        
      }
    }
  }
}
#endif
