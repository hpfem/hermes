#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "h1.h"
#include "../function/forms.h"
#include <algorithm>

namespace WeakFormsNeutronics
{
  namespace MaterialProperties
  { 
    typedef double rank0;
    typedef std::vector<double> rank1;
    typedef std::vector<std::vector<double > > rank2;
    typedef std::vector<std::vector<std::vector<double > > > rank3;
    
    typedef std::map<std::string, rank0> MaterialPropertyMap0;
    typedef std::map<std::string, rank1> MaterialPropertyMap1;
    typedef std::map<std::string, rank2> MaterialPropertyMap2;
    typedef std::map<std::string, rank3> MaterialPropertyMap3;
    
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
    static const char* E_NONMATCHING_PROPERTIES =
      "All properties must be defined for a single given number of materials.";
    static const char* E_INVALID_SIZE =
      "Material property defined for an unexpected number of groups.";
    static const char* E_INVALID_GROUP_INDEX =
      "Attempted to access property data in an out-of-range group.";
    static const char* E_SG_SIGMA_R = 
      "Group-reduction cross-section (Sigma_r) is not defined for one-group (i.e. monoenergetic) problems.";
      
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
                typedef typename NDArrayType::value_type dim_type;                      \
                typename NDArrayType::const_iterator dim_iterator_x = x.begin();        \
                typename NDArrayType::const_iterator dim_iterator_y = y.begin();        \
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
    
    struct ensure_trivial { 
      void operator() (MaterialPropertyMap1::value_type x) { 
        MaterialPropertyMap1::mapped_type::iterator it;
        for (it = x.second.begin(); it != x.second.end(); ++it) 
          if (fabs(*it) > 1e-14)
            error(E_INVALID_COMBINATION);
      }
    };
    
    struct ensure_size { 
      ensure_size(unsigned int nrows, unsigned int ncols = -1) 
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
        
        MaterialPropertyMaps(unsigned int G, std::set<std::string> mat_list = std::set<std::string>()) 
          : materials_list(mat_list), G(G)  { };
                     
        virtual void validate()
        {                                              
          if (nu.empty() && !nuSigma_f.empty() && !Sigma_f.empty())
            nu = map_divide<rank1>(nuSigma_f, Sigma_f);
          else if (nuSigma_f.empty() && !nu.empty() && !Sigma_f.empty())
            nuSigma_f = map_multiply<rank1>(nu, Sigma_f);
          else if (Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
            Sigma_f = map_divide<rank1>(nuSigma_f, nu);
          else if (!Sigma_f.empty() && !nuSigma_f.empty() && !nu.empty())
          {
            MaterialPropertyMap1 diff = map_subtract<rank1>(nuSigma_f, map_multiply<rank1>(nu, Sigma_f));            
            std::for_each(diff.begin(), diff.end(), ensure_trivial());
          }
          else
            error(E_INSUFFICIENT_DATA);
          
          if ((nu.size() != Sigma_f.size()) || (nu.size() != chi.size()))
            error(E_NONMATCHING_PROPERTIES);
          
          std::for_each(nu.begin(), nu.end(), ensure_size(G));
          std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size(G));
          std::for_each(chi.begin(), chi.end(), ensure_size(G));
          
          if (Sigma_a.size() > 0)
          {
            // Warn if \Sigma_a < \Sigma_f for any region (this indicates an unphysical situation, since
            // by definition \Sigma_a = \Sigma_f + \Sigma_c + \Sigma_{n,p} + other possible reactions
            // leading to neutron removal).
            MaterialPropertyMap1::const_iterator ita = Sigma_a.begin();
            MaterialPropertyMap1::const_iterator itf = Sigma_f.begin();
            for ( ; ita != Sigma_a.end(); ++ita, ++itf)
              if (*ita < *itf)
                warning(W_SA_LT_SF);
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
    using namespace MaterialProperties;
    
    typedef std::map<int,rank0> iMarkerPropertyMap0;
    typedef std::map<int,rank1> iMarkerPropertyMap1;
    typedef std::map<int,rank2> iMarkerPropertyMap2;
    typedef std::map<int,rank3> iMarkerPropertyMap3;
    
    struct ensure_size_at_least { 
      ensure_size_at_least(unsigned int nrows, unsigned int ncols = -1) 
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
        
    class GenericMultigroupWeakForm : public WeakForm
    {
      protected:
                
        iMarkerPropertyMap1 Sigma_f,
                            nu,
                            chi;
      
      public:
        
        GenericMultigroupWeakForm(unsigned int G, const MaterialPropertyMaps& matprop) : WeakForm(G)
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
        }
    };
              
    namespace Diffusion
    {
      class DiffusionMaterialPropertyMaps : public MaterialPropertyMaps
      {
        private:
                                      
          MaterialPropertyMap1 D;
          MaterialPropertyMap1 Sigma_r;
          MaterialPropertyMap2 Sigma_s;
          MaterialPropertyMap1 src;
          
          MaterialPropertyMap1 Sigma_t;
                                    
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
            
            if ((D.size() != Sigma_r.size()) || (D.size() != Sigma_s.size()) || (D.size() != src.size()))
              error(E_NONMATCHING_PROPERTIES);
            
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
      
      class GenericMultigroupDiffusionWeakForm : public GenericMultigroupWeakForm
      {
        protected:
                    
          iMarkerPropertyMap1 D, Sigma_r;
          iMarkerPropertyMap2 Sigma_s;
        
        public:
          
          GenericMultigroupDiffusionWeakForm(unsigned int G, const DiffusionMaterialPropertyMaps& matprop) 
            : GenericMultigroupWeakForm(G, matprop)
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
            
            MaterialPropertyMap2 str_mpm2 = matprop.get_Sigma_s();
            MaterialPropertyMap2::const_iterator str_mpm2_iter = str_mpm2.begin();
            for ( ; str_mpm2_iter != str_mpm2.end(); ++str_mpm2_iter)
            {
              int int_marker = element_markers_conversion->get_internal_marker(str_mpm2_iter->first);
              Sigma_s[int_marker] = str_mpm2_iter->second;
            }
          }
      };
      
      namespace MatrixForms 
      {  
        class DefaultDiffusionReaction : public WeakForm::MatrixFormVol
        {
          public:
            
            DefaultDiffusionReaction(unsigned int g, 
                                    const iMarkerPropertyMap1& D, const iMarkerPropertyMap1& Sigma_r,
                                    GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(g, g, HERMES_ANY, HERMES_SYM),
                g(g), D(D), Sigma_r(Sigma_r), geom_type(geom_type)
            {
              if (D.size() != Sigma_r.size()) 
                error(E_NONMATCHING_PROPERTIES);
              
              std::for_each(D.begin(), D.end(), ensure_size_at_least(g));
              std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size_at_least(g));
            }
            DefaultDiffusionReaction(unsigned int g, std::string area, 
                                    const iMarkerPropertyMap1& D, const iMarkerPropertyMap1& Sigma_r,
                                    GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(g, g, area, HERMES_SYM),
                g(g), D(D), Sigma_r(Sigma_r), geom_type(geom_type)
            { 
              if (D.size() != Sigma_r.size())
                error(E_NONMATCHING_PROPERTIES);
              
              std::for_each(D.begin(), D.end(), ensure_size_at_least(g));
              std::for_each(Sigma_r.begin(), Sigma_r.end(), ensure_size_at_least(g));
            }
            
            template<typename Real, typename Scalar>
            Scalar matrix_form( int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext  ) const 
            {
              Scalar result = 0;
              
              // Constant properties within the current active element. Note that prop[e->elem_marker]
              // cannot be used since 'prop' is a constant std::map for which operator[] is undefined. 
              rank1 D_elem = D.find(e->elem_marker)->second;
              rank1 Sigma_r_elem = Sigma_r.find(e->elem_marker)->second;
              
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
              return new DefaultDiffusionReaction(*this);
            }

          private:
            
            unsigned int g;
            iMarkerPropertyMap1 D, Sigma_r;
            GeomType geom_type;
        };
                
        class DefaultFissionYield : public WeakForm::MatrixFormVol
        {
          public:
            
            DefaultFissionYield(unsigned int gto, unsigned int gfrom, 
                                const iMarkerPropertyMap1& chi,
                                const iMarkerPropertyMap1& Sigma_f,
                                const iMarkerPropertyMap1& nu,
                                GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(gto, gfrom), 
                gto(gto), gfrom(gfrom), chi(chi), Sigma_f(Sigma_f), nu(nu), geom_type(geom_type)
            {
              if (chi.size() != Sigma_f.size() || chi.size() != nu.size())
                error(E_NONMATCHING_PROPERTIES);
              
              std::for_each(chi.begin(), chi.end(), ensure_size_at_least(gto));
              std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size_at_least(gfrom));
              std::for_each(nu.begin(), nu.end(), ensure_size_at_least(gfrom));
            }
            
            DefaultFissionYield(unsigned int gto, unsigned int gfrom, std::string area,
                                const iMarkerPropertyMap1& chi,
                                const iMarkerPropertyMap1& Sigma_f,
                                const iMarkerPropertyMap1& nu,
                                GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(gto, gfrom, area),
                gto(gto), gfrom(gfrom), chi(chi), Sigma_f(Sigma_f), nu(nu), geom_type(geom_type)
            { 
              if (chi.size() != Sigma_f.size() || chi.size() != nu.size())
                error(E_NONMATCHING_PROPERTIES);
              
              std::for_each(chi.begin(), chi.end(), ensure_size_at_least(gto));
              std::for_each(Sigma_f.begin(), Sigma_f.end(), ensure_size_at_least(gfrom));
              std::for_each(nu.begin(), nu.end(), ensure_size_at_least(gfrom));
            }
            
            virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                  Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const 
            {
              scalar result = 0;
              if (geom_type == HERMES_PLANAR) result = int_u_v<double, scalar>(n, wt, u, v);
              else 
              {
                if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<double, scalar>(n, wt, u, v, e);
                else result = int_x_u_v<double, scalar>(n, wt, u, v, e);
              }
              
              // Constant properties within the current active element. Note that prop[e->elem_marker]
              // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
              rank1 nu_elem = nu.find(e->elem_marker)->second;
              rank1 Sigma_f_elem = Sigma_f.find(e->elem_marker)->second;
              rank1 chi_elem = chi.find(e->elem_marker)->second;
              
              return result * chi_elem[gto] * nu_elem[gfrom] * Sigma_f_elem[gfrom];
            }
            
            // This is to make the form usable in rk_time_step().
            virtual WeakForm::MatrixFormVol* clone() {
              return new DefaultFissionYield(*this);
            }
            
          private:
            
            unsigned int gto, gfrom;
            iMarkerPropertyMap1 chi, Sigma_f, nu;
            GeomType geom_type;
        };
            
        class DefaultScattering : public WeakForm::MatrixFormVol
        {
          public:
            
            DefaultScattering(unsigned int gto, unsigned int gfrom, 
                              const iMarkerPropertyMap2& Sigma_s,
                              GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(gto, gfrom), 
                gto(gto), gfrom(gfrom), Sigma_s(Sigma_s), geom_type(geom_type)
            {
              std::for_each(Sigma_s.begin(), Sigma_s.end(), ensure_size_at_least(gto, gfrom));
            }
            
            DefaultScattering(unsigned int gto, unsigned int gfrom, std::string area,
                              const iMarkerPropertyMap2& Sigma_s,
                              GeomType geom_type = HERMES_PLANAR)
              : WeakForm::MatrixFormVol(gto, gfrom, area),
                gto(gto), gfrom(gfrom), Sigma_s(Sigma_s), geom_type(geom_type)
            { 
              std::for_each(Sigma_s.begin(), Sigma_s.end(), ensure_size_at_least(gto, gfrom));
            }
            
            virtual scalar value( int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                  Func<double> *v, Geom<double> *e, ExtData<scalar> *ext  ) const 
            {
              scalar result = 0;
              if (geom_type == HERMES_PLANAR) result = int_u_v<double, scalar>(n, wt, u, v);
              else 
              {
                if (geom_type == HERMES_AXISYM_X) result = int_y_u_v<double, scalar>(n, wt, u, v, e);
                else result = int_x_u_v<double, scalar>(n, wt, u, v, e);
              }
              
              // Constant properties within the current active element. Note that prop[e->elem_marker]
              // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
              rank2 Sigma_s_elem = Sigma_s.find(e->elem_marker)->second;
              return result * (-Sigma_s_elem[gto][gfrom]);
            }
            
            // This is to make the form usable in rk_time_step().
            virtual WeakForm::MatrixFormVol* clone() {
              return new DefaultScattering(*this);
            }
            
          private:
            
            unsigned int gto, gfrom;
            iMarkerPropertyMap2 Sigma_s;
            GeomType geom_type;
        };
      };
      
      namespace VectorForms 
      {        
        class DefaultFissionYieldIterative : public WeakForm::VectorFormVol
        {
          public:
            
            DefaultFissionYieldIterative(unsigned int g, 
                                         const iMarkerPropertyMap1& chi,
                                         const iMarkerPropertyMap1& Sigma_f,
                                         const iMarkerPropertyMap1& nu,
                                         int keff = 1.0,
                                         GeomType geom_type = HERMES_PLANAR)
              : WeakForm::VectorFormVol(g),
                g(g), keff(keff), chi(chi), Sigma_f(Sigma_f), nu(nu), geom_type(geom_type)
            {
              if (chi.size() != Sigma_f.size() || chi.size() != nu.size())
                error(E_NONMATCHING_PROPERTIES);
              std::for_each(chi.begin(), chi.end(), ensure_size_at_least(g));
            }
            
            DefaultFissionYieldIterative(unsigned int g, std::string area,
                                         const iMarkerPropertyMap1& chi,
                                         const iMarkerPropertyMap1& Sigma_f,
                                         const iMarkerPropertyMap1& nu,
                                         int keff = 1.0,
                                         GeomType geom_type = HERMES_PLANAR)
              : WeakForm::VectorFormVol(g, area),
                g(g), keff(keff), chi(chi), Sigma_f(Sigma_f), nu(nu), geom_type(geom_type)
            {
              if (chi.size() != Sigma_f.size() || chi.size() != nu.size())
                error(E_NONMATCHING_PROPERTIES);
              std::for_each(chi.begin(), chi.end(), ensure_size_at_least(g));
            }
            
            template<typename Real, typename Scalar>
            Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
            { 
              // Constant properties within the current active element. Note that prop[e->elem_marker]
              // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
              rank1 nu_elem = nu.find(e->elem_marker)->second;
              rank1 Sigma_f_elem = Sigma_f.find(e->elem_marker)->second;
              rank1 chi_elem = chi.find(e->elem_marker)->second;
              
              if ((unsigned)ext->nf != nu_elem.size() || (unsigned)ext->nf != Sigma_f.size())
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
                                Geom<double> *e, ExtData<scalar> *ext) const {
              return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
            }

            virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
              return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
            }

            // This is to make the form usable in rk_time_step().
            virtual WeakForm::VectorFormVol* clone() {
              return new DefaultFissionYieldIterative(*this);
            }
            
            void update_keff(double new_keff) { keff = new_keff; }
            
          private:
            
            unsigned int g;
            double keff;
            iMarkerPropertyMap1 chi, Sigma_f, nu;
            GeomType geom_type;
        };
        
        class DefaultExternalSource : public WeakForm::VectorFormVol
        {
          public:
            
            DefaultExternalSource(unsigned int g, 
                                  const iMarkerPropertyMap1& src,
                                  GeomType geom_type = HERMES_PLANAR)
              : WeakForm::VectorFormVol(g), g(g), src(src), geom_type(geom_type) 
            { 
              std::for_each(src.begin(), src.end(), ensure_size_at_least(g));
            }
            
            DefaultExternalSource(unsigned int g, std::string area,
                                  const iMarkerPropertyMap1& src,
                                  GeomType geom_type = HERMES_PLANAR)
              : WeakForm::VectorFormVol(g, area), g(g), src(src), geom_type(geom_type)
            { 
              std::for_each(src.begin(), src.end(), ensure_size_at_least(g));
            }
            
            virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const 
            {
              // Constant properties within the current active element. Note that prop[e->elem_marker]
              // cannot be used since 'prop' is a constant std::map for which operator[] is undefined.
              rank1 src_elem = src.find(e->elem_marker)->second;
              
              if (geom_type == HERMES_PLANAR) return src_elem[g] * int_v<double>(n, wt, v);
              else 
              {
                if (geom_type == HERMES_AXISYM_X) 
                  return src_elem[g] * int_y_v<double>(n, wt, v, e);
                else 
                  return src_elem[g] * int_x_v<double>(n, wt, v, e);
              }
            }
            
            // This is to make the form usable in rk_time_step().
            virtual WeakForm::VectorFormVol* clone() {
              return new DefaultExternalSource(*this);
            }
          
          private:
            
            unsigned int g;
            iMarkerPropertyMap1 src;
            GeomType geom_type;        
        };
      }
          
      class DefaultWeakFormFixedSource : public GenericMultigroupDiffusionWeakForm
      { 
        DefaultWeakFormFixedSource(unsigned int G, const DiffusionMaterialPropertyMaps& matprop) 
          : GenericMultigroupDiffusionWeakForm(G, matprop)
        {
          
        }
      };
      
      class DefaultWeakFormSourceIteration : public GenericMultigroupDiffusionWeakForm
      {
      };
      
    };
  }
}
#endif