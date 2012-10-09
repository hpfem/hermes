#ifndef __H2D_NEUTRONICS_WEAK_FORMS_H
#define __H2D_NEUTRONICS_WEAK_FORMS_H

#include "weakforms_h1.h"
#include "../forms.h"
#include "../function/filter.h"

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
          template<typename Scalar>
          class HERMES_API DefaultWeakFormFixedSource : public WeakForm<Scalar>
          {
          public:
            DefaultWeakFormFixedSource( Hermes::vector<std::string> regions,
              Hermes::vector<double> D_map,
              Hermes::vector<double> Sigma_a_map,
              Hermes::vector<double> Q_map );
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

            typedef std::vector<bool > bool1;
            typedef std::vector<std::vector<bool > > bool2;
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
                  if(fabs(*it) > 1e-14)
                    throw Hermes::Exceptions::Exception(E_INVALID_COMBINATION);
              }
            };

            struct ensure_size {
              ensure_size(unsigned int nrows, unsigned int ncols = 0)
                : nrows(nrows), ncols(ncols) {};

              void operator() (MaterialPropertyMap1::value_type x) {
                if(x.second.size() != nrows)
                  throw Hermes::Exceptions::Exception(E_INVALID_SIZE);
              }

              void operator() (MaterialPropertyMap2::value_type x) {
                if(x.second.size() != nrows)
                  throw Hermes::Exceptions::Exception(E_INVALID_SIZE);

                MaterialPropertyMap2::mapped_type::iterator it;
                for (it = x.second.begin(); it != x.second.end(); ++it)
                  if(it->size() != ncols)
                    throw Hermes::Exceptions::Exception(E_INVALID_SIZE);
              }

            private:
              unsigned int nrows, ncols;
            };
          }

          namespace Common
          {
            using namespace Definitions;
            using namespace Messages;

            class HERMES_API NDArrayMapOp
            {
              //
              // NOTE: Could be perhaps combined with the classes material_property_map and MultiArray below
              // and moved to hermes_common as a general way of handling maps with multidimensional mapped types.
              //

              template <typename NDArrayType>
              static rank0 divide(rank0 x, rank0 y) {
                if(x == 0 && y == 0)
                  return 0.0;
                else if(y == 0)
                {
                  throw Hermes::Exceptions::Exception(E_INF_VALUE);
                  return -1.0;
                }
                else
                  return x/y;
              }

              template <typename NDArrayType>
              static rank0 multiply(rank0 x, rank0 y) {
                return x*y;
              }

              template <typename NDArrayType>
              static rank0 add(rank0 x, rank0 y) {
                return x + y;
              }

              template <typename NDArrayType>
              static rank0 subtract(rank0 x, rank0 y) {
                rank0 ret = x - y;
                return ret;
              }

            public:

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

            class HERMES_API MaterialPropertyMaps : public Hermes::Mixins::Loggable
            {
            protected:

              MaterialPropertyMap1 Sigma_f;
              MaterialPropertyMap1 nu;
              MaterialPropertyMap1 chi;

              MaterialPropertyMap1 Sigma_a;
              MaterialPropertyMap1 nuSigma_f;

              std::set<std::string> materials_list;
              unsigned int G;

              bool1 fission_multigroup_structure;

              void extend_to_multigroup(const MaterialPropertyMap0& mrsg_map, MaterialPropertyMap1 *mrmg_map);

              void extend_to_multiregion(const rank1& srmg_array, MaterialPropertyMap1 *mrmg_map);

              void extend_to_multiregion_multigroup(const rank0& srsg_value, MaterialPropertyMap1 *mrmg_map);

              void fill_with(double c, MaterialPropertyMap1 *mrmg_map);

              MaterialPropertyMaps(unsigned int G, std::set<std::string> mat_list = std::set<std::string>())
                : materials_list(mat_list), G(G)  { };

              virtual void validate();

            public:

              void set_nu(const MaterialPropertyMap1& nu) {
                this->nu = nu;
              }

              void set_nu(const MaterialPropertyMap0& nu) {
                extend_to_multigroup(nu, &this->nu);
              }

              void set_nu(const rank1& nu) {
                extend_to_multiregion(nu, &this->nu);
              }

              void set_nu(const rank0& nu) {
                extend_to_multiregion_multigroup(nu, &this->nu);
              }

              void set_chi(const MaterialPropertyMap1& chi) {
                this->chi = chi;
              }

              void set_chi(const rank1& chi) {
                extend_to_multiregion(chi, &this->chi);
              }

              void set_fission_multigroup_structure(const bool1& chi_nnz)  {
                this->fission_multigroup_structure = chi_nnz;
              }

              void set_Sigma_a(const MaterialPropertyMap1& Sa) {
                this->Sigma_a = Sa;
              }

              void set_Sigma_f(const MaterialPropertyMap1& Sf) {
                this->Sigma_f = Sf;
              }

              void set_nuSigma_f(const MaterialPropertyMap1 nSf) {
                this->nuSigma_f = nSf;
              }

              const MaterialPropertyMap1& get_Sigma_f() const {
                return this->Sigma_f;
              }
              const MaterialPropertyMap1& get_nu() const {
                return this->nu;
              }
              const MaterialPropertyMap1& get_chi() const {
                return this->chi;
              }
              const bool1& get_fission_multigroup_structure() const {
                return this->fission_multigroup_structure;
              }

              const rank1& get_Sigma_f(std::string material) const;
              const rank1& get_nu(std::string material) const;
              const rank1& get_chi(std::string material) const;

              unsigned int get_G() const { return G; }

              friend std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
            };
          }

          namespace Diffusion
          {
            using namespace Definitions;
            using namespace Messages;

            class HERMES_API MaterialPropertyMaps : public Common::MaterialPropertyMaps
            {
            protected:

              MaterialPropertyMap1 D;
              MaterialPropertyMap1 Sigma_r;
              MaterialPropertyMap2 Sigma_s;
              MaterialPropertyMap1 src;

              MaterialPropertyMap1 Sigma_t;

              bool2 scattering_multigroup_structure;

            public:

              MaterialPropertyMaps(unsigned int G, std::set<std::string> mat_list = std::set<std::string>())
                : Common::MaterialPropertyMaps(G, mat_list) { };

              MaterialPropertyMap1 extract_map2_diagonals(const MaterialPropertyMap2& map2);

              MaterialPropertyMap1 sum_map2_columns(const MaterialPropertyMap2& map2);

              MaterialPropertyMap2 create_map2_by_diagonals(const MaterialPropertyMap1& diags);

              void fill_with(double c, MaterialPropertyMap2 *mrmg_map);

              // We always need to supply chi, nu, Sigma_f, Sigma_r, Sigma_s and D to our neutronics weak forms.
              // These parameters are often defined in terms of the other ones, or not specified at all and assumed
              // to be zero for a particular simplified situation. This method, together with its complement in the
              // parent class, uses the most typical definitions to build the six-parameter set from the given input.
              // It also checks whether the user did not enter nonsensical values. However, values entered by the
              // user may sometimes not satisfy the common relations, as some empirical corrections may have been
              // already included in them.
              virtual void validate();

              void set_src(const MaterialPropertyMap1& src) {
                this->src = src;
              }

              void set_src(const MaterialPropertyMap0& src) {
                extend_to_multigroup(src, &this->src);
              }

              void set_src(const rank1& src) {
                extend_to_multiregion(src, &this->src);
              }

              void set_src(const double& src) {
                extend_to_multiregion_multigroup(src, &this->src);
              }

              void set_D(const MaterialPropertyMap1& D) {
                this->D = D;
              }

              void set_Sigma_r(const MaterialPropertyMap1& Sr) {
                this->Sigma_r = Sr;
              }

              void set_Sigma_t(const MaterialPropertyMap1& St) {
                this->Sigma_t = St;
              }

              void set_Sigma_s(const MaterialPropertyMap2& Ss) {
                this->Sigma_s = Ss;
              }

              void set_scattering_multigroup_structure(const bool2& Ss_nnz) {
                this->scattering_multigroup_structure = Ss_nnz;
              }

              const MaterialPropertyMap2& get_Sigma_s() const {
                return this->Sigma_s;
              }
              const MaterialPropertyMap1& get_Sigma_r() const {
                return this->Sigma_r;
              }
              const MaterialPropertyMap1& get_D() const {
                return this->D;
              }
              const MaterialPropertyMap1& get_src() const {
                return this->src;
              }
              const bool2& get_scattering_multigroup_structure() const {
                return this->scattering_multigroup_structure;
              }

              const rank2& get_Sigma_s(std::string material) const;
              const rank1& get_Sigma_r(std::string material) const;
              const rank1& get_D(std::string material) const;
              const rank1& get_src(std::string material) const;

              friend HERMES_API std::ostream & operator<< (std::ostream& os, const MaterialPropertyMaps& matprop);
            };
          }

          template <typename NDArrayType>
          class material_property_map
          {
          private:
            std::map<std::string, NDArrayType> m_map;
          public:
            material_property_map(const std::string& key, const NDArrayType& val) {
              m_map[key] = val;
            }

            material_property_map<NDArrayType>& operator()(const std::string& key, const NDArrayType& val) {
              m_map[key] = val;
              return *this;
            }

            operator std::map<std::string, NDArrayType>() {
              return m_map;
            }
          };

          template <typename NDArrayType>
          class MultiArray
          {
          private:
            Hermes::vector<NDArrayType> m_data;
          public:
            MultiArray(const NDArrayType& val) {
              m_data.push_back(val);
            }

            MultiArray<NDArrayType>& operator()(const NDArrayType& val) {
              m_data.push_back(val);
              return *this;
            }

            operator Hermes::vector<NDArrayType>() {
              return m_data;
            }
          };

          namespace Definitions
          {
            typedef MultiArray<rank0> row;
            typedef MultiArray<rank1> mat;
            typedef MultiArray<bool> bool_row;
            typedef MultiArray< Hermes::vector<bool> > bool_mat;
          }
        }

        namespace ElementaryForms
        {
          namespace Diffusion
          {
            using namespace MaterialProperties::Diffusion;

            class HERMES_API GenericForm
            {
            protected:
              const MaterialPropertyMaps& matprop;
              GeomType geom_type;

              GenericForm(const MaterialPropertyMaps& matprop,
                GeomType geom_type = HERMES_PLANAR)
                : matprop(matprop), geom_type(geom_type)
              {};

              GenericForm(const MaterialPropertyMaps& matprop, Mesh* mesh,
                GeomType geom_type = HERMES_PLANAR)
                : matprop(matprop), geom_type(geom_type), mesh(mesh)
              {};

              Mesh* mesh;
            };

            struct HERMES_API VacuumBoundaryCondition
            {
              // TODO: General albedo boundary condition.
              template<typename Scalar>
              class HERMES_API Jacobian : public MatrixFormSurf<Scalar>
              {
              public:
                Jacobian(unsigned int g, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormSurf<Scalar>(g, g),
                  g(g), geom_type(geom_type)
                {};

                Jacobian(unsigned int g, std::string area, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormSurf<Scalar>(g, g),
                  g(g), geom_type(geom_type)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns matrix_form(int n, double *wt, Func<ScalarTestFns> *u_ext[], Func<Real> *u,
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                  Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const {
                    return matrix_form<double, Scalar>(n, wt, u_ext, u, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
                  Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const {
                    return matrix_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual MatrixFormSurf<Scalar>* clone() {
                  return new Jacobian(*this);
                }

              private:
                unsigned int g;
                GeomType geom_type;
              };

              template<typename Scalar>
              class HERMES_API Residual : public VectorFormSurf<Scalar>
              {
              public:
                Residual(unsigned int g, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormSurf<Scalar>(g),
                  g(g), geom_type(geom_type)
                {};

                Residual(unsigned int g, std::string area, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormSurf<Scalar>(g),
                  g(g), geom_type(geom_type)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
                  Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const {
                    return vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
                  Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormSurf<Scalar>* clone() {
                  return new Residual(*this);
                }

              private:
                unsigned int g;
                GeomType geom_type;
              };
            };

            struct HERMES_API DiffusionReaction
            {
              template<typename Scalar>
              class Jacobian : public MatrixFormVol<Scalar>, protected GenericForm
              {
              public:
                Jacobian(unsigned int g,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(g, g),
                  GenericForm(matprop, geom_type),
                  g(g)
                {};

                Jacobian(unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(g, g),
                  GenericForm(matprop, geom_type),
                  g(g)
                {
                  this->set_area(area);
                };

                Jacobian(unsigned int g,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(g, g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {};

                Jacobian(unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(g, g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns matrix_form( int n, double *wt, Func<ScalarTestFns> *u_ext[], Func<Real> *u,
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext  ) const;

                virtual Scalar value( int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                  Func<double> *v, Geom<double> *e, Func<Scalar> **ext  ) const {
                    return  matrix_form<double, Scalar> (n, wt, u_ext, u, v, e, ext);
                }

                virtual Hermes::Ord ord( int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
                  Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Hermes::Ord> **ext  ) const {
                    return  matrix_form<Hermes::Ord, Hermes::Ord> (n, wt, u_ext, u, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual MatrixFormVol<Scalar>* clone() {
                  return new Jacobian(*this);
                }

              private:

                unsigned int g;
              };

              template<typename Scalar>
              class Residual : public VectorFormVol<Scalar>, protected GenericForm
              {
              public:

                Residual(unsigned int g,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, geom_type),
                  g(g)
                {};

                Residual(unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, geom_type),
                  g(g)
                {
                  this->set_area(area);
                };

                Residual(unsigned int g,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {};

                Residual(unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                  Geom<double> *e, Func<Scalar> **ext) const  {
                    return vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, Func<Ord> **ext) const  {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormVol<Scalar>* clone() {
                  return new Residual(*this);
                }

              private:

                unsigned int g;
              };
            };

            struct HERMES_API FissionYield
            {
              template<typename Scalar>
              class HERMES_API Jacobian : public MatrixFormVol<Scalar>, protected GenericForm
              {
              public:

                Jacobian( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {
                  this->set_area(area);
                };

                Jacobian( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, Mesh* mesh, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh* mesh, GeomType geom_type = HERMES_PLANAR)
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns matrix_form( int n, double *wt, Func<ScalarTestFns> *u_ext[], Func<Real> *u,
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext  ) const;

                virtual Scalar value( int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                  Func<double> *v, Geom<double> *e, Func<Scalar> **ext  ) const {
                    return  -1.0 * matrix_form<double, Scalar> (n, wt, u_ext, u, v, e, ext);
                }

                virtual Hermes::Ord ord( int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
                  Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Hermes::Ord> **ext  ) const {
                    return  matrix_form<Hermes::Ord, Hermes::Ord> (n, wt, u_ext, u, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual MatrixFormVol<Scalar>* clone() {
                  return new Jacobian(*this);
                }

              private:

                unsigned int gto, gfrom;
              };

              template<typename Scalar>
              class HERMES_API OuterIterationForm : public VectorFormVol<Scalar>, protected GenericForm
              {
              public:

                OuterIterationForm( unsigned int g,
                  const MaterialPropertyMaps& matprop,
                  Hermes::vector<MeshFunction<Scalar>*>& iterates,
                  double keff = 1.0,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, geom_type),
                  g(g), keff(keff)
                {
                  this->wf->set_ext(iterates);
                  if(g >= iterates.size())
                    throw Hermes::Exceptions::Exception(E_INVALID_GROUP_INDEX);
                }

                OuterIterationForm( unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop,
                  Hermes::vector<MeshFunction<Scalar>*>& iterates,
                  double keff = 1.0,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, geom_type),
                  g(g), keff(keff)
                {
                  this->set_area(area);
                  this->wf->set_ext(iterates);
                  if(g >= iterates.size())
                    throw Hermes::Exceptions::Exception(E_INVALID_GROUP_INDEX);
                }

                OuterIterationForm( unsigned int g,
                  const MaterialPropertyMaps& matprop, Mesh *mesh,
                  Hermes::vector<MeshFunction<Scalar>*>& iterates,
                  double keff = 1.0,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g), keff(keff)
                {
                  this->wf->set_ext(iterates);
                  if(g >= iterates.size())
                    throw Hermes::Exceptions::Exception(E_INVALID_GROUP_INDEX);
                }

                OuterIterationForm( unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh,
                  Hermes::vector<MeshFunction<Scalar>*>& iterates,
                  double keff = 1.0,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g), keff(keff)
                {
                  this->set_area(area);
                  this->wf->set_ext(iterates);
                  if(g >= iterates.size())
                    throw Hermes::Exceptions::Exception(E_INVALID_GROUP_INDEX);
                }

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                  Geom<double> *e, Func<Scalar> **ext) const {
                    return -1.0 * vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, Func<Ord> **ext) const  {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormVol<Scalar>* clone() {
                  return new OuterIterationForm(*this);
                }

                void update_keff(double new_keff) { keff = new_keff; }

              private:

                unsigned int g;
                double keff;
              };

              template<typename Scalar>
              class HERMES_API Residual : public VectorFormVol<Scalar>, protected GenericForm
              {
              public:
                Residual( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {
                  this->set_area(area);
                };

                Residual( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {
                  this->set_area(area);
                };

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                  Geom<double> *e, Func<Scalar> **ext) const {
                    return -1.0 * vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, Func<Ord> **ext) const {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormVol<Scalar>* clone() {
                  return new Residual(*this);
                }

              private:

                unsigned int gto, gfrom;
              };
            };

            struct HERMES_API Scattering
            {
              template<typename Scalar>
              class Jacobian : public MatrixFormVol<Scalar>, protected GenericForm
              {
              public:

                Jacobian( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom, area),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Jacobian( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Jacobian( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR )
                  : MatrixFormVol<Scalar>(gto, gfrom, area),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns matrix_form( int n, double *wt, Func<ScalarTestFns> *u_ext[], Func<Real> *u,
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext  ) const;

                virtual Scalar value( int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
                  Func<double> *v, Geom<double> *e, Func<Scalar> **ext  ) const {
                    return  -1.0 * matrix_form<double, Scalar> (n, wt, u_ext, u, v, e, ext);
                }

                virtual Hermes::Ord ord( int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
                  Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Hermes::Ord> **ext  ) const {
                    return  matrix_form<Hermes::Ord, Hermes::Ord> (n, wt, u_ext, u, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual MatrixFormVol<Scalar>* clone() {
                  return new Jacobian(*this);
                }

              private:

                unsigned int gto, gfrom;
              };

              template<typename Scalar>
              class Residual : public VectorFormVol<Scalar>, protected GenericForm
              {
              public:
                Residual( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto, area),
                  GenericForm(matprop, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Residual( unsigned int gto, unsigned int gfrom,
                  const MaterialPropertyMaps& matprop, Mesh *mesh,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                Residual( unsigned int gto, unsigned int gfrom, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh,
                  GeomType geom_type = HERMES_PLANAR )
                  : VectorFormVol<Scalar>(gto, area),
                  GenericForm(matprop, mesh, geom_type),
                  gto(gto), gfrom(gfrom)
                {};

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                  Geom<double> *e, Func<Scalar> **ext) const {
                    return -1.0 * vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, Func<Ord> **ext) const {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormVol<Scalar>* clone() {
                  return new Residual(*this);
                }

              private:

                unsigned int gto, gfrom;
              };
            };

            struct HERMES_API ExternalSources
            {
              template<typename Scalar>
              class LinearForm : public VectorFormVol<Scalar>, protected GenericForm
              {
              public:

                LinearForm( unsigned int g,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, geom_type),
                  g(g)
                {};

                LinearForm( unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g, area),
                  GenericForm(matprop, geom_type),
                  g(g)
                {};

                LinearForm( unsigned int g,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {};

                LinearForm( unsigned int g, std::string area,
                  const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type = HERMES_PLANAR)
                  : VectorFormVol<Scalar>(g, area),
                  GenericForm(matprop, mesh, geom_type),
                  g(g)
                {};

                template<typename Real, typename ScalarTestFns>
                ScalarTestFns vector_form(int n, double *wt, Func<ScalarTestFns> *u_ext[],
                  Func<Real> *v, Geom<Real> *e, Func<ScalarTestFns> **ext) const;

                virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
                  Geom<double> *e, Func<Scalar> **ext) const {
                    return -1.0 * vector_form<double, Scalar>(n, wt, u_ext, v, e, ext);
                }

                virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, Func<Ord> **ext) const {
                    return vector_form<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
                }

                // This is to make the form usable in rk_time_step_newton().
                virtual VectorFormVol<Scalar>* clone() {
                  return new LinearForm(*this);
                }

              private:

                unsigned int g;
              };
            };
          }
        }

        namespace CompleteWeakForms
        {
          namespace Diffusion
          {
            using namespace MaterialProperties;
            using namespace ElementaryForms::Diffusion;

            template<typename Scalar>
            class HERMES_API DefaultWeakFormFixedSource : public WeakForm<Scalar>
            {
            protected:
              void lhs_init(unsigned int G, const MaterialPropertyMaps& matprop, Mesh *mesh, GeomType geom_type);

            public:
              DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, Mesh *mesh,
                GeomType geom_type = HERMES_PLANAR);

              DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, Mesh *mesh,
                Hermes2DFunction<Scalar>*f_src,
                std::string src_area = HERMES_ANY,
                GeomType geom_type = HERMES_PLANAR);

              DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, Mesh *mesh,
                Hermes2DFunction<Scalar>*f_src,
                Hermes::vector<std::string> src_areas,
                GeomType geom_type = HERMES_PLANAR);

              DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, Mesh *mesh,
                const Hermes::vector<Hermes2DFunction<Scalar>*>& f_src,
                std::string src_area = HERMES_ANY,
                GeomType geom_type = HERMES_PLANAR);

              DefaultWeakFormFixedSource(const MaterialPropertyMaps& matprop, Mesh *mesh,
                const Hermes::vector<Hermes2DFunction<Scalar>*>& f_src,
                Hermes::vector<std::string> src_areas,
                GeomType geom_type = HERMES_PLANAR);
            };

            template<typename Scalar>
            class HERMES_API DefaultWeakFormSourceIteration : public WeakForm<Scalar>
            {
            protected:
              Hermes::vector<FissionYield::OuterIterationForm<Scalar>*> keff_iteration_forms;

            public:
              DefaultWeakFormSourceIteration( const MaterialPropertyMaps& matprop, Mesh *mesh,
                Hermes::vector<MeshFunction<Scalar>*>& iterates,
                double initial_keff_guess,
                GeomType geom_type = HERMES_PLANAR );

              void update_keff(double new_keff);
              /// \todo This is needed by 4-group adapt, however it must have been removed, so I provided this dummy method to
              /// get over it.
              double get_keff() { return 0.0; };
            };
          }
        }

        namespace SupportClasses
        {
          using MaterialProperties::Common::MaterialPropertyMaps;
          using namespace MaterialProperties::Definitions;

          class HERMES_API SourceFilter : public SimpleFilter<double>
          {
            public:
            SourceFilter(Hermes::vector<MeshFunction<double>*> solutions, const MaterialPropertyMaps* matprop,
                         const std::string& source_area)
              : SimpleFilter<double>(solutions, Hermes::vector<int>())
            {
              nu = matprop->get_nu().at(source_area);
              Sigma_f = matprop->get_Sigma_f().at(source_area);
            }
            SourceFilter(Hermes::vector<Solution<double>*> solutions, const MaterialPropertyMaps* matprop,
                const std::string& source_area)
            : SimpleFilter<double>(solutions, Hermes::vector<int>())
            {
              nu = matprop->get_nu().at(source_area);
              Sigma_f = matprop->get_Sigma_f().at(source_area);
            }

          private:
            rank1 nu;
            rank1 Sigma_f;

            void filter_fn(int n, Hermes::vector<double*> values, double* result);
          };
        }
      }
    }
  }
}
#endif