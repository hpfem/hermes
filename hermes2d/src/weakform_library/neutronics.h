#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "h1.h"
#include <algorithm>

/* Default weak form for neutron diffusion equations
   with Dirichlet and/or zero Neumann BC.

   Nonzero Neumann or Newton boundary conditions can be enabled
   by creating a descendant and adding surface forms to it.
*/

namespace WeakFormsNeutronics
{
  namespace MaterialProperties
  { 
    typedef double rank0;
    typedef std::vector<double> rank1;
    typedef std::vector<std::vector<double > > rank2;
    typedef std::vector<std::vector<std::vector<double > > > rank3;
    
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
    static const char* E_INVALID_COMBINATION =
      "Invalid combination of entered material properties.";
          
    template <typename NDArrayType>
    rank0 divide(rank0 x, rank0 y) 
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
    rank0 multiply(rank0 x, rank0 y) 
    {
      return x*y;
    }
    
    template <typename NDArrayType>
    rank0 add(rank0 x, rank0 y) 
    {
      return x + y;
    }
    
    template <typename NDArrayType>
    rank0 subtract(rank0 x, rank0 y) 
    {
      rank0 ret = x - y;
      if(ret < 0)
        warning(W_NEG_VALUE);
      return ret;
    }
    
    #define for_each_element_in_dimension \
                typedef typename NDArrayType::value_type dim_type;                        \
                typename NDArrayType::const_iterator dim_iterator_x = x.begin();          \
                typename NDArrayType::const_iterator dim_iterator_y = y.begin();          \
                for ( ; dim_iterator_x != x.end(); ++dim_iterator_x, ++dim_iterator_y )                    
    
    template <typename NDArrayType>
    NDArrayType divide(const NDArrayType& x, const NDArrayType& y) 
    { 
      NDArrayType res; res.reserve(x.size());
            
      for_each_element_in_dimension
        res.push_back( divide<dim_type>(*dim_iterator_x, *dim_iterator_x) );
      
      return res;
    }
    
    template <typename NDArrayType>
    NDArrayType multiply(const NDArrayType& x, const NDArrayType& y) 
    { 
      NDArrayType res; res.reserve(x.size());
      
      for_each_element_in_dimension
        res.push_back( multiply<dim_type>(*dim_iterator_x, *dim_iterator_x) );
      
      return res;
    }
    
    template <typename NDArrayType>
    NDArrayType add(const NDArrayType& x, const NDArrayType& y) 
    { 
      NDArrayType res; res.reserve(x.size());
      
      for_each_element_in_dimension
        res.push_back( add<dim_type>(*dim_iterator_x, *dim_iterator_x) );
      
      return res;
    }
    
    template <typename NDArrayType>
    NDArrayType subtract(const NDArrayType& x, const NDArrayType& y) 
    { 
      NDArrayType res; res.reserve(x.size());
      
      for_each_element_in_dimension
        res.push_back( subtract<dim_type>(*dim_iterator_x, *dim_iterator_x) );
      
      return res;
    }
    
    
    #undef for_each_element_in_dimension
    
    #define for_each_element_in_map \
              typename std::map<std::string, T>::iterator iterator_ret = ret.begin();   \
              typename std::map<std::string, T>::const_iterator iterator_x = x.begin(); \
              typename std::map<std::string, T>::const_iterator iterator_y = y.begin(); \
              for ( ; iterator_x != x.end(); ++iterator_x, ++iterator_y, ++iterator_ret )                    
    
    template <typename T>
    std::map<std::string,  T> map_divide(const std::map<std::string, T>& x, 
                                         const std::map<std::string, T>& y)
    {
      std::map<std::string, T> ret = x;
      
      for_each_element_in_map
        iterator_ret->second = divide<T>(iterator_x->second, iterator_y->second);
      
      return ret;
    }
    
    template <typename T>
    std::map<std::string,  T> map_multiply(const std::map<std::string, T>& x, 
                                           const std::map<std::string, T>& y)
    {
      std::map<std::string, T> ret = x;
      
      for_each_element_in_map
        iterator_ret->second = multiply<T>(iterator_x->second, iterator_y->second);
      
      return ret;
    }
    
    template <typename T>
    std::map<std::string,  T> map_add(const std::map<std::string, T>& x, 
                                      const std::map<std::string, T>& y)
    {
      std::map<std::string, T> ret = x;
      
      for_each_element_in_map
        iterator_ret->second = add<T>(iterator_x->second, iterator_y->second);
      
      return ret;
    }
    
    template <typename T>
    std::map<std::string,  T> map_subtract(const std::map<std::string, T>& x, 
                                           const std::map<std::string, T>& y)
    {
      std::map<std::string, T> ret = x;
      
      for_each_element_in_map
        iterator_ret->second = subtract<T>(iterator_x->second, iterator_y->second);
      
      return ret;
    }
    
    #undef for_each_element_in_map
    
    class MaterialPropertyMaps
    {
      protected:
        
        typedef std::map<std::string, rank0> MatPropMap0;
        typedef std::map<std::string, rank1> MatPropMap1;
        
        MatPropMap1 Sigma_f;
        MatPropMap1 nu;
        MatPropMap1 chi;
        
        MatPropMap1 Sigma_a;
        MatPropMap1 nuSigma_f;
        
        std::set<std::string> materials_list;
        int G;
              
        void extend_to_multigroup(const MatPropMap0& mrsg_map, MatPropMap1 *mrmg_map)
        {
          if (G == 1)
            warning(W_MG_EXTENSION);
          
          MatPropMap0::const_iterator it;
          for (it = mrsg_map.begin(); it != mrsg_map.end(); ++it)
            mrmg_map->at(it->first).assign(G, it->second);
        
        }
        
        void extend_to_multiregion(const rank1& srmg_array, MatPropMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            mrmg_map->at(*it) = srmg_array;
        }
        
        void extend_to_multiregion_multigroup(const rank0& srsg_value, MatPropMap1 *mrmg_map)
        {
          if (materials_list.empty())
            error(E_MR_EXTENSION);
          if (G == 1)
            warning(W_MG_EXTENSION);
          
          std::set<std::string>::const_iterator it;
          for (it = materials_list.begin(); it != materials_list.end(); ++it)
            mrmg_map->at(*it).assign(G, srsg_value);
        }
        
        MaterialPropertyMaps(int G, std::set<std::string> mat_list = std::set<std::string>()) 
          : materials_list(mat_list), G(G)  { };
          
        struct ensure_trivial_rank1 { 
          void operator() (MatPropMap1::value_type x) { 
            MatPropMap1::mapped_type::iterator it;
            for (it = x.second.begin(); it != x.second.end(); ++it) 
              if (fabs(*it) > 1e-14)
                error(E_INVALID_COMBINATION);
          }
        };
            
        virtual void validate()
        {
          // Warn if \Sigma_a < \Sigma_f for any region (this indicates an unphysical situation, since
          // by definition \Sigma_a = \Sigma_f + \Sigma_c + \Sigma_{n,p} + other possible reactions
          // leading to neutron removal).
          MatPropMap1::const_iterator ita = Sigma_a.begin();
          MatPropMap1::const_iterator itf = Sigma_f.begin();
          for ( ; ita != Sigma_a.end(); ++ita, ++itf)
            if (*ita < *itf)
              warning(W_SA_LT_SF);
                        
          if (nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
            nu = map_divide<rank1>(nuSigma_f, Sigma_f);
          else if (nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
            nuSigma_f = map_multiply<rank1>(nu, Sigma_f);
          else if (Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
            Sigma_f = map_divide<rank1>(nuSigma_f, nu);
          else if (!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
          {
            MatPropMap1 diff = map_subtract<rank1>(nuSigma_f, map_multiply<rank1>(nu, Sigma_f));            
            std::for_each(diff.begin(), diff.end(), ensure_trivial_rank1());
          }
          else
            error(E_INSUFFICIENT_DATA);
        }
        
      public:
        
        void set_nu(const MatPropMap1& nu)
        {
          this->nu = nu;
        }
        
        void set_nu(const MatPropMap0& nu)
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
              
        void set_chi(const MatPropMap1& nu)
        {
          this->chi = chi;
        }
        
        void set_chi(const rank1& chi)
        {
          extend_to_multiregion(chi, &this->chi);
        }
        
        void set_Sigma_a(const MatPropMap1& Sa)
        {
          this->Sigma_a = Sa;
        }
        
        void set_Sigma_f(const MatPropMap1& Sf)
        {
          this->Sigma_f = Sf;
        }
        
        void set_nuSigma_f(const MatPropMap1 nSf)
        {
          this->nuSigma_f = nSf;
        }  
    };
    
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
    
    typedef MultiArray<rank0> grow;
    typedef MultiArray<rank1> gmat;
  };

  namespace Diffusion 
  {    
    namespace Monoenergetic
    {
      /* Simple monoenergetic neutron diffusion, with the following weak formulation within each
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
      */
      class DefaultWeakFormFixedSource : public WeakForm
      {        
        public:
          DefaultWeakFormFixedSource( Hermes::vector<std::string> regions, 
                                      Hermes::vector<double> D_map, 
                                      Hermes::vector<double> Sigma_a_map, 
                                      Hermes::vector<double> Q_map ) : WeakForm(1) 
          {
            using namespace WeakFormsH1::VolumetricMatrixForms;
            using namespace WeakFormsH1::VolumetricVectorForms;
            
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
    };
    
    namespace Multigroup
    {   
      using namespace MaterialProperties;
      
      static const char* E_SG_SIGMA_R = 
        "Group-reduction cross-section (Sigma_r) is not defined for one-group (i.e. monoenergetic) problems.";
      
      class DiffusionMaterialPropertyMaps : public MaterialPropertyMaps
      {
        private:
          
          typedef std::map<std::string, rank2> MatPropMap2;
                            
          MatPropMap1 D;
          MatPropMap1 Sigma_r;
          MatPropMap2 Sigma_s;
          MatPropMap1 src;
          
          MatPropMap1 Sigma_t;
                                    
        public:
          
          DiffusionMaterialPropertyMaps(int G, std::set<std::string> mat_list = std::set<std::string>()) 
            : MaterialPropertyMaps(G, mat_list) { };
                  
          virtual void validate()
          {
            MaterialPropertyMaps::validate();
            
            if (G == 1) // Monoenergetic case.
            {    
            }
            else  // Multigroup case.
            {
            }
          }
          
          void set_src(const MatPropMap1& src)
          {
            this->src = src;
          }
          
          void set_src(const MatPropMap0& src)
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
          
          void set_D(const MatPropMap1& D)
          {
            this->D = D;
          }
          
          void set_Sigma_r(const MatPropMap1& Sr)
          {
            this->Sigma_r = Sr;
          }
          
          void set_Sigma_t(const MatPropMap1& St)
          {
            this->Sigma_t = St;
          }
          
          void set_Sigma_s(const MatPropMap2& Ss)
          {
            this->Sigma_s = Ss;
          }      
      };
      
      namespace MatrixForms 
      {  
        class DefaultDiffusionReaction : public WeakForm::MatrixFormVol
        {
        };
        
        class DefaultFissionYield : public WeakForm::MatrixFormVol
        {
        };
            
        class DefaultScattering : public WeakForm::MatrixFormVol
        {
        };
      };
      
      namespace VectorForms 
      {
        class DefaultFissionYieldIterative : public WeakForm::VectorFormVol
        {
        };
      }
          
      class DefaultWeakFormFixedSource : public WeakForm
      {
        DefaultWeakFormFixedSource(int G, MaterialPropertyMaps matprop) : WeakForm(G)
        {
        }
      };
      
      class DefaultWeakFormSourceIteration : public WeakForm
      {
      };
      
    };  
  }
}
#endif