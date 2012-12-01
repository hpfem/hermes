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
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ExactSolution<Scalar>::ExactSolution(const Mesh* mesh) : Solution<Scalar>(mesh)
    {
      this->sln_type = HERMES_EXACT;
      this->num_dofs = -1;
      this->exact_multiplicator = 1.0;
    }

    template<typename Scalar>
    MeshFunction<Scalar>* ExactSolution<Scalar>::clone() const
    {
      throw Hermes::Exceptions::Exception("Solution<Scalar>::clone() must be overridden in the case of exact solutions.");
      return NULL;
    }

    template<typename Scalar>
    ExactSolutionScalar<Scalar>::ExactSolutionScalar(const Mesh* mesh) : ExactSolution<Scalar>(mesh)
    {
      this->num_components = 1;
    }

    template<typename Scalar>
    unsigned int ExactSolutionScalar<Scalar>::get_dimension() const
    {
      return 1;
    }

    template<typename Scalar>
    ExactSolutionVector<Scalar>::ExactSolutionVector(const Mesh* mesh) : ExactSolution<Scalar>(mesh)
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

        xmlsolution.exactConstantXReal() = this->constant;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
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

        xmlsolution.exactConstantXReal() = this->constant.real();
        xmlsolution.exactConstantXComplex() = this->constant.imag();

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<typename Scalar>
    ConstantSolution<Scalar>::ConstantSolution(const Mesh* mesh, Scalar constant) : ExactSolutionScalar<Scalar>(mesh), constant(constant) {};

    template<typename Scalar>
    Scalar ConstantSolution<Scalar>::value (double x, double y) const {
      return constant;
    };

    template<typename Scalar>
    MeshFunction<Scalar>* ConstantSolution<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      ConstantSolution<Scalar>* sln = new ConstantSolution<Scalar>(this->mesh, this->constant);
      return sln;
    }

    template<typename Scalar>
    void ConstantSolution<Scalar>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const {
      dx = 0;
      dy = 0;
    };

    template<typename Scalar>
    Ord ConstantSolution<Scalar>::ord(Ord x, Ord y) const {
      return Ord(0);
    }

    template<>
    void ZeroSolution<double>::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(1, 0, 0, 1, 0);

        xmlsolution.exactConstantXReal() = 0;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void ZeroSolution<std::complex<double> >::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(1, 0, 0, 1, 1);

        xmlsolution.exactConstantXReal() = 0;
        xmlsolution.exactConstantXComplex() = 0;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<typename Scalar>
    ZeroSolution<Scalar>::ZeroSolution(const Mesh* mesh) : ExactSolutionScalar<Scalar>(mesh) {};

    template<typename Scalar>
    Scalar ZeroSolution<Scalar>::value (double x, double y) const {
      return 0.0;
    };

    template<typename Scalar>
    MeshFunction<Scalar>* ZeroSolution<Scalar>::clone() const
    {
      if(this->sln_type == HERMES_SLN)
        return Solution<Scalar>::clone();
      ZeroSolution<Scalar>* sln = new ZeroSolution<Scalar>(this->mesh);
      return sln;
    }

    template<typename Scalar>
    void ZeroSolution<Scalar>::derivatives (double x, double y, Scalar& dx, Scalar& dy) const {
      dx = 0;
      dy = 0;
    };

    template<typename Scalar>
    Ord ZeroSolution<Scalar>::ord(Ord x, Ord y) const {
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

        xmlsolution.exactConstantXReal() = this->constantX;
        xmlsolution.exactConstantYReal() = this->constantY;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
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

        xmlsolution.exactConstantXReal() = this->constantX.real();
        xmlsolution.exactConstantXComplex() = this->constantX.imag();
        xmlsolution.exactConstantYReal() = this->constantY.real();
        xmlsolution.exactConstantYComplex() = this->constantY.imag();

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<typename Scalar>
    ConstantSolutionVector<Scalar>::ConstantSolutionVector(const Mesh* mesh, Scalar constantX, Scalar constantY) : ExactSolutionVector<Scalar>(mesh), constantX(constantX), constantY(constantY) {};

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
    Ord ConstantSolutionVector<Scalar>::ord(Ord x, Ord y) const {
      return Ord(0);
    }

    template<>
    void ZeroSolutionVector<double>::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<double>::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(2, 0, 0, 1, 0);

        xmlsolution.exactConstantXReal() = 0;
        xmlsolution.exactConstantYReal() = 0;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<>
    void ZeroSolutionVector<std::complex<double> >::save(const char* filename) const
    {
      if(this->sln_type == HERMES_SLN)
      {
        Solution<std::complex<double> >::save(filename);
        return;
      }
      try
      {
        XMLSolution::solution xmlsolution(2, 0, 0, 1, 1);

        xmlsolution.exactConstantXReal() = 0;
        xmlsolution.exactConstantXComplex() = 0;
        xmlsolution.exactConstantYReal() = 0;
        xmlsolution.exactConstantYComplex() = 0;

        std::string solution_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
        solution_schema_location.append("/solution_h2d_xml.xsd");
        ::xml_schema::namespace_info namespace_info_solution("XMLSolution", solution_schema_location);

        ::xml_schema::namespace_infomap namespace_info_map;
        namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("solution", namespace_info_solution));

        std::ofstream out(filename);
        XMLSolution::solution_(out, xmlsolution, namespace_info_map);
        out.close();
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SolutionSaveFailureException(e.what());
      }
    }

    template<typename Scalar>
    ZeroSolutionVector<Scalar>::ZeroSolutionVector(const Mesh* mesh) : ExactSolutionVector<Scalar>(mesh) {};

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
    Ord ZeroSolutionVector<Scalar>::ord(Ord x, Ord y) const {
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

    template HERMES_API class ExactSolutionScalar<double>;
    template HERMES_API class ExactSolutionScalar<std::complex<double> >;
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
  }
}
