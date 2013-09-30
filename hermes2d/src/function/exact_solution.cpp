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

#include "solution_h2d_xml.h"
#include "exact_solution.h"
#include "api2d.h"
#include "../weakform_library/weakforms_h1.h"
#include "space_h1.h"
#include "../solver/linear_solver.h"

#ifdef WITH_BSON
#include "bson.h"
#endif

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ExactSolution<Scalar>::ExactSolution(MeshSharedPtr mesh) : Solution<Scalar>(mesh)
    {
      this->sln_type = HERMES_EXACT;
      this->num_dofs = -1;
      this->exact_multiplicator = 1.0;
    }

    template<typename Scalar>
    void ExactSolution<Scalar>::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<Scalar>::save(filename);
        return;
      }

      throw Exceptions::Exception("Arbitrary exact solution can not be saved to disk. Only constant one can. Project to a space to get a saveable solution.");
    }

#ifdef WITH_BSON
    template<typename Scalar>
    void ExactSolution<Scalar>::save_bson(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<Scalar>::save_bson(filename);
        return;
      }

      throw Exceptions::Exception("Arbitrary exact solution can not be saved to disk. Only constant one can. Project to a space to get a saveable solution.");
    }
#endif

    template<typename Scalar>
    MeshFunction<Scalar>* ExactSolution<Scalar>::clone() const
    {
      throw Hermes::Exceptions::Exception("Solution<Scalar>::clone() must be overridden in the case of exact solutions.");
      return NULL;
    }

    template<typename Scalar>
    ExactSolutionScalar<Scalar>::ExactSolutionScalar(MeshSharedPtr mesh) : ExactSolution<Scalar>(mesh)
    {
      this->num_components = 1;
    }

    template<typename Scalar>
    unsigned int ExactSolutionScalar<Scalar>::get_dimension() const
    {
      return 1;
    }

    template<typename Scalar, typename ValueType>
    ExactSolutionConstantArray<Scalar, ValueType>::ExactSolutionConstantArray(MeshSharedPtr mesh, ValueType* valueArray, bool deleteArray) : ExactSolutionScalar<Scalar>(mesh), valueArray(valueArray), deleteArray(deleteArray)
    {
    }

    template<typename Scalar, typename ValueType>
    ExactSolutionConstantArray<Scalar, ValueType>::~ExactSolutionConstantArray()
    {
      if(this->deleteArray)
        delete [] this->valueArray;
    }

    template<typename Scalar, typename ValueType>
    MeshFunction<Scalar>* ExactSolutionConstantArray<Scalar, ValueType>::clone() const
    {
      return new ExactSolutionConstantArray<Scalar, ValueType>(this->mesh, this->valueArray);
    }

    template<typename Scalar, typename ValueType>
    Scalar ExactSolutionConstantArray<Scalar, ValueType>::value (double x, double y) const
    {
      return this->valueArray[this->element->id];
    }

    template<typename Scalar, typename ValueType>
    Ord ExactSolutionConstantArray<Scalar, ValueType>::ord(double x, double y) const {
      return Ord(0);
    }

    template<typename Scalar, typename ValueType>
    void ExactSolutionConstantArray<Scalar, ValueType>::setArray(ValueType* valueArray)
    {
      valueArray = valueArray;
    }

    template<typename Scalar, typename ValueType>
    void ExactSolutionConstantArray<Scalar, ValueType>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const {
      dx = 0;
      dy = 0;
    };



    template<typename Scalar>
    ExactSolutionVector<Scalar>::ExactSolutionVector(MeshSharedPtr mesh) : ExactSolution<Scalar>(mesh)
    {
      this->num_components = 2;
    }

    template<typename Scalar>
    unsigned int ExactSolutionVector<Scalar>::get_dimension() const
    {
      return 2;
    }

    template<>
    void ConstantSolution<double>::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(1, 0, 0, 1, 0);

        xmlsolution.exactCXR() = this->constant;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;
        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);

        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void ConstantSolution<std::complex<double> >::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(1, 0, 0, 1, 1);

        xmlsolution.exactCXR() = this->constant.real();
        xmlsolution.exactCXC() = this->constant.imag();

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;
        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

#ifdef WITH_BSON
    template<>
    void ConstantSolution<double>::save_bson(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save_bson(filename);
        return;
      }

      // bson
      bson bw;
      bson_init(&bw);

      bson_append_bool(&bw, "exact", true);
      bson_append_bool(&bw, "complex", false);

      bson_append_start_array(&bw, "values");
      bson_append_double(&bw, "c", this->constant);
      bson_append_double(&bw, "c", 0.);
      bson_append_double(&bw, "c", 0.);
      bson_append_double(&bw, "c", 0.);
      bson_append_finish_array(&bw);

      bson_append_int(&bw, "components_count", this->num_components);

      bson_finish(&bw);

      FILE *fpw;
      fpw = fopen(filename, "wb");
      const char *dataw = (const char *) bson_data(&bw);
      fwrite(dataw, bson_size(&bw), 1, fpw);
      fclose(fpw);

      bson_destroy(&bw);
    }

    template<>
    void ConstantSolution<std::complex<double> >::save_bson(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save_bson(filename);
        return;
      }

      // bson
      bson bw;
      bson_init(&bw);

      bson_append_bool(&bw, "exact", true);
      bson_append_bool(&bw, "complex", true);

      bson_append_start_array(&bw, "values");
      bson_append_double(&bw, "c", this->constant.real());
      bson_append_double(&bw, "c", 0.);
      bson_append_double(&bw, "c", this->constant.imag());
      bson_append_double(&bw, "c", 0.);
      bson_append_finish_array(&bw);

      bson_append_int(&bw, "components_count", this->num_components);

      bson_finish(&bw);

      FILE *fpw;
      fpw = fopen(filename, "wb");
      const char *dataw = (const char *) bson_data(&bw);
      fwrite(dataw, bson_size(&bw), 1, fpw);
      fclose(fpw);

      bson_destroy(&bw);
    }
#endif

    template<typename Scalar>
    ConstantSolution<Scalar>::ConstantSolution(MeshSharedPtr mesh, Scalar constant) : ExactSolutionScalar<Scalar>(mesh), constant(constant)
    {
      this->order = 0;
    };

    template<typename Scalar>
    Scalar ConstantSolution<Scalar>::value (double x, double y) const {
      return constant;
    };

    template<typename Scalar>
    MeshFunction<Scalar>* ConstantSolution<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      return new ConstantSolution<Scalar>(this->mesh, this->constant);
    }

    template<typename Scalar>
    void ConstantSolution<Scalar>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const {
      dx = 0;
      dy = 0;
    };

    template<typename Scalar>
    Ord ConstantSolution<Scalar>::ord(double x, double y) const {
      return Ord(0);
    }

    template<typename Scalar>
    ZeroSolution<Scalar>::ZeroSolution(MeshSharedPtr mesh) : ExactSolutionScalar<Scalar>(mesh)
    {
      this->order = 0;
    };

    template<typename Scalar>
    Scalar ZeroSolution<Scalar>::value (double x, double y) const {
      return 0.0;
    };

    template<typename Scalar>
    MeshFunction<Scalar>* ZeroSolution<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      return new ZeroSolution<Scalar>(this->mesh);
    }

    template<typename Scalar>
    void ZeroSolution<Scalar>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const {
      dx = 0;
      dy = 0;
    };

    template<typename Scalar>
    Ord ZeroSolution<Scalar>::ord(double x, double y) const {
      return Ord(0);
    }


    template<>
    void ConstantSolutionVector<double>::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(2, 0, 0, 1, 0);

        xmlsolution.exactCXR() = this->constantX;
        xmlsolution.exactCYR() = this->constantY;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;
        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void ConstantSolutionVector<std::complex<double> >::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(2, 0, 0, 1, 1);

        xmlsolution.exactCXR() = this->constantX.real();
        xmlsolution.exactCXC() = this->constantX.imag();
        xmlsolution.exactCYR() = this->constantY.real();
        xmlsolution.exactCYC() = this->constantY.imag();

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);

        ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;

        XMLSolution::solution_(out, xmlsolution, namespace_info_map, "UTF-8", parsing_flags);

        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<typename Scalar>
    ConstantSolutionVector<Scalar>::ConstantSolutionVector(MeshSharedPtr mesh, Scalar constantX, Scalar constantY) : ExactSolutionVector<Scalar>(mesh), constantX(constantX), constantY(constantY)
    {
      this->order = 0;
    };

    template<typename Scalar>
    MeshFunction<Scalar>* ConstantSolutionVector<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      ConstantSolutionVector<Scalar>* sln = new ConstantSolutionVector<Scalar>(this->mesh, this->constantX, this->constantY);
      return sln;
    }

    template<typename Scalar>
    Scalar2<Scalar> ConstantSolutionVector<Scalar>::value (double x, double y) const {
      return Scalar2<Scalar>(constantX, constantY);
    };

    template<typename Scalar>
    void ConstantSolutionVector<Scalar>::derivatives (double x, double y, Scalar2<Scalar>& dx, Scalar2<Scalar>& dy) const {
      dx = Scalar2<Scalar>(Scalar(0.0), Scalar(0.0));
      dy = Scalar2<Scalar>(Scalar(0.0), Scalar(0.0));
    };

    template<typename Scalar>
    Ord ConstantSolutionVector<Scalar>::ord(double x, double y) const {
      return Ord(0);
    }

#ifdef WITH_BSON
    template<>
    void ConstantSolutionVector<double>::save_bson(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save_bson(filename);
        return;
      }

      // bson
      bson bw;
      bson_init(&bw);

      bson_append_bool(&bw, "exact", true);
      bson_append_bool(&bw, "complex", false);

      bson_append_start_array(&bw, "values");
      bson_append_double(&bw, "c", this->constantX);
      bson_append_double(&bw, "c", this->constantY);
      bson_append_double(&bw, "c", 0.);
      bson_append_double(&bw, "c", 0.);
      bson_append_finish_array(&bw);

      bson_append_int(&bw, "components_count", this->num_components);

      bson_finish(&bw);

      FILE *fpw;
      fpw = fopen(filename, "wb");
      const char *dataw = (const char *) bson_data(&bw);
      fwrite(dataw, bson_size(&bw), 1, fpw);
      fclose(fpw);

      bson_destroy(&bw);
    }

    template<>
    void ConstantSolutionVector<std::complex<double> >::save_bson(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save_bson(filename);
        return;
      }

      // bson
      bson bw;
      bson_init(&bw);

      bson_append_bool(&bw, "exact", true);
      bson_append_bool(&bw, "complex", true);

      bson_append_start_array(&bw, "values");
      bson_append_double(&bw, "c", this->constantX.real());
      bson_append_double(&bw, "c", this->constantY.real());
      bson_append_double(&bw, "c", this->constantX.imag());
      bson_append_double(&bw, "c", this->constantY.imag());
      bson_append_finish_array(&bw);

      bson_append_int(&bw, "components_count", this->num_components);

      bson_finish(&bw);

      FILE *fpw;
      fpw = fopen(filename, "wb");
      const char *dataw = (const char *) bson_data(&bw);
      fwrite(dataw, bson_size(&bw), 1, fpw);
      fclose(fpw);

      bson_destroy(&bw);
    }
#endif

    template<typename Scalar>
    ZeroSolutionVector<Scalar>::ZeroSolutionVector(MeshSharedPtr mesh) : ExactSolutionVector<Scalar>(mesh)
    {
      this->order = 0;
    };

    template<typename Scalar>
    Scalar2<Scalar> ZeroSolutionVector<Scalar>::value (double x, double y) const {
      return Scalar2<Scalar>(0.0, 0.0);
    };

    template<typename Scalar>
    void ZeroSolutionVector<Scalar>::derivatives (double x, double y, Scalar2<Scalar>& dx, Scalar2<Scalar>& dy) const {
      dx = Scalar2<Scalar>(0.0, 0.0);
      dy = Scalar2<Scalar>(0.0, 0.0);
    };

    template<typename Scalar>
    Ord ZeroSolutionVector<Scalar>::ord(double x, double y) const {
      return Ord(0);
    }

    template<typename Scalar>
    MeshFunction<Scalar>* ZeroSolutionVector<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      ZeroSolutionVector<Scalar>* sln = new ZeroSolutionVector<Scalar>(this->mesh);
      return sln;
    }

    ExactSolutionEggShell::ExactSolutionEggShell(MeshSharedPtr mesh, int polynomialOrder) : ExactSolutionScalar<double>(mesh)
    {
      Hermes2D::WeakFormsH1::DefaultWeakFormLaplaceLinear<double> wf;
      SpaceSharedPtr<double> space(new H1SpaceEggShell(mesh, polynomialOrder));
      Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
      MeshFunctionSharedPtr<double> sln(new Solution<double>());
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space, sln);
      this->copy(sln.get());
    }

    double ExactSolutionEggShell::value (double x, double y) const
    {
      throw Exceptions::Exception("ExactSolutionEggShell::value should never be called.");
      return 0.;
    }

    void ExactSolutionEggShell::derivatives (double x, double y, double& dx, double& dy) const
    {
      throw Exceptions::Exception("ExactSolutionEggShell::derivatives should never be called.");
    }

    Hermes::Ord ExactSolutionEggShell::ord(double x, double y) const
    {
      throw Exceptions::Exception("ExactSolutionEggShell::ord should never be called.");
      return Hermes::Ord(0);
    }

    MeshFunction<double>* ExactSolutionEggShell::clone() const
    {
      Solution<double> * sln = new Solution<double>;
      sln->copy(this);
      return sln;
    }

    template<typename Scalar>
    UExtFunction<Scalar>::UExtFunction(MeshSharedPtr mesh) : MeshFunction<Scalar>(mesh)
    {
      this->num_components = 1;
      this->set_quad_2d(&g_quad_2d_std);
    }

    template HERMES_API class ExactSolutionScalar<double>;
    template HERMES_API class ExactSolutionScalar<std::complex<double> >;

    template HERMES_API class ExactSolutionConstantArray<double, double>;
    template HERMES_API class ExactSolutionConstantArray<double, int>;
    template HERMES_API class ExactSolutionConstantArray<double, unsigned int>;
    template HERMES_API class ExactSolutionConstantArray<double, bool>;
    template HERMES_API class ExactSolutionConstantArray<std::complex<double>, std::complex<double> >;

    template HERMES_API class ExactSolutionVector<double>;
    template HERMES_API class ExactSolutionVector<std::complex<double> >;
    template HERMES_API class ConstantSolution<double>;
    template HERMES_API class ConstantSolution<std::complex<double> >;
    template HERMES_API class ConstantSolutionVector<double>;
    template HERMES_API class ConstantSolutionVector<std::complex<double> >;
    template HERMES_API class ZeroSolution<double>;
    template HERMES_API class ZeroSolution<std::complex<double> >;
    template HERMES_API class ZeroSolutionVector<double>;
    template HERMES_API class ZeroSolutionVector<std::complex<double> >;

    template HERMES_API class UExtFunction<double>;
    template HERMES_API class UExtFunction<std::complex<double> >;
  }
}
