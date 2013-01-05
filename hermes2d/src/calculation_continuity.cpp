// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include "calculation_continuity.h"
#include "mesh_reader_h2d_xml.h"
#include "space_h1.h"
#include "space_hdiv.h"
#include "space_hcurl.h"
#include "space_l2.h"
#include <fstream>

namespace Hermes
{
  namespace Hermes2D
  {
    CalculationContinuityException::CalculationContinuityException() : Exception()
    {
    }

    CalculationContinuityException::CalculationContinuityException(exceptionEntityType type, const char * reason) : Exception()
    {
      this->init(type, reason);
    }

    void CalculationContinuityException::init(exceptionEntityType type, const char * reason)
    {
      char * msg =  new char[34 + strlen(reason)];
      char * typeMsg = new char[15];
      switch(type)
      {
      case meshes:
        sprintf(typeMsg, "meshes");
        break;
      case spaces:
        sprintf(typeMsg, "spaces");
        break;
      case solutions:
        sprintf(typeMsg, "solutions");
        break;
      case time_steps:
        sprintf(typeMsg, "time steps");
        break;
      case error:
        sprintf(typeMsg, "error info");
        break;
      case general:
        sprintf(typeMsg, "general");
        break;
      }

      sprintf(msg, "Exception in CalculationContinuity (%s): \"%s\"", typeMsg, reason);
      message = msg;
      delete [] typeMsg;
    }

    IOCalculationContinuityException::IOCalculationContinuityException(exceptionEntityType type, inputOutput inputOutput, const char * filename) : CalculationContinuityException()
    {
      char * msg =  new char[34 + strlen(filename)];
      char * typeMsg = new char[7];
      switch(inputOutput)
      {
      case input:
        sprintf(typeMsg, "input");
        break;
      case output:
        sprintf(typeMsg, "output");
        break;
      }
      sprintf(msg, "I/O exception: %s, filename: \"%s\"", typeMsg, filename);
      this->init(type, msg);
    }

    IOCalculationContinuityException::IOCalculationContinuityException(exceptionEntityType type, inputOutput inputOutput, const char * filename, const char * reason) : CalculationContinuityException()
    {
      char * msg =  new char[34 + strlen(filename)];
      char * typeMsg = new char[7];
      switch(inputOutput)
      {
      case input:
        sprintf(typeMsg, "input");
        break;
      case output:
        sprintf(typeMsg, "output");
        break;
      }
      sprintf(msg, "I/O exception: %s, filename: \"%s\", reason: %s", typeMsg, filename, reason);
      this->init(type, msg);
    }

    template<typename Scalar>
    CalculationContinuity<Scalar>::CalculationContinuity(IdentificationMethod identification_method) : last_record(NULL), record_available(false), identification_method(identification_method), num(0)
    {
      double last_time;
      unsigned int last_number;
      std::stringstream ss;
      switch(identification_method)
      {
      case timeAndNumber:
        ss << "timeAndNumber.h2d";
        break;
      case onlyTime:
        ss << "onlyTime.h2d";
        break;
      case onlyNumber:
        ss << "onlyNumber.h2d";
        break;
      }

      std::ifstream ifile(ss.str().c_str());
      if(ifile)
      {
        while(!ifile.eof())
        {
          switch(identification_method)
          {
          case timeAndNumber:
            ifile >> this->num >> last_time >> last_number;
            record_available = true;
            break;
          case onlyTime:
            ifile >> this->num >> last_time;
            record_available = true;
            break;
          case onlyNumber:
            ifile >> this->num >> last_number;
            record_available = true;
            break;
          }
        }
        ifile.close();

        CalculationContinuity<Scalar>::Record* record;
        switch(identification_method)
        {
        case timeAndNumber:
          record = new CalculationContinuity<Scalar>::Record(last_time, last_number);
          this->records.insert(std::pair<std::pair<double, unsigned int>, CalculationContinuity<Scalar>::Record*>(std::pair<double, unsigned int>(last_time, last_number), record));
          this->last_record = record;
          break;
        case onlyTime:
          record = new CalculationContinuity<Scalar>::Record(last_time);
          this->time_records.insert(std::pair<double, CalculationContinuity<Scalar>::Record*>(last_time, record));
          this->last_record = record;
          break;
        case onlyNumber:
          record = new CalculationContinuity<Scalar>::Record(last_number);
          this->numbered_records.insert(std::pair<unsigned int, CalculationContinuity<Scalar>::Record*>(last_number, record));
          this->last_record = record;
          break;
        }
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(double time, unsigned int number, Mesh* mesh, Space<Scalar>* space, Solution<Scalar>* sln, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("timeAndNumber.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << time << ' ' << number << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "timeAndNumber.h2d");

      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(time, number);
      record->save_mesh(mesh);
      if(space != NULL)
        record->save_space(space);
      if(sln != NULL)
        record->save_solution(sln);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);

      this->records.insert(std::pair<std::pair<double, unsigned int>, CalculationContinuity<Scalar>::Record*>(std::pair<double, unsigned int>(time, number), record));
      this->last_record = record;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(double time, unsigned int number, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> slns, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("timeAndNumber.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << time << ' ' << number << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "timeAndNumber.h2d");

      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(time, number);
      record->save_meshes(meshes);
      if(spaces != Hermes::vector<Space<Scalar>*>())
        record->save_spaces(spaces);
      if(slns != Hermes::vector<Solution<Scalar>*>())
        record->save_solutions(slns);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);

      this->records.insert(std::pair<std::pair<double, unsigned int>, CalculationContinuity<Scalar>::Record*>(std::pair<double, unsigned int>(time, number), record));
      this->last_record = record;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(double time, Mesh* mesh, Space<Scalar>* space, Solution<Scalar>* sln, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("onlyTime.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << time << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "onlyTime.h2d");

      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(time);
      record->save_mesh(mesh);
      if(space != NULL)
        record->save_space(space);
      if(sln != NULL)
        record->save_solution(sln);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);

      this->time_records.insert(std::pair<double, CalculationContinuity<Scalar>::Record*>(time, record));
      this->last_record = record;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(double time, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> slns, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("onlyTime.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << time << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "onlyTime.h2d");
      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(time);
      record->save_meshes(meshes);
      if(spaces != Hermes::vector<Space<Scalar>*>())
        record->save_spaces(spaces);
      if(slns != Hermes::vector<Solution<Scalar>*>())
        record->save_solutions(slns);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);
      this->time_records.insert(std::pair<double, CalculationContinuity<Scalar>::Record*>(time, record));
      this->last_record = record;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(unsigned int number, Mesh* mesh, Space<Scalar>* space, Solution<Scalar>* sln, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("onlyNumber.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << number << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "onlyNumber.h2d");

      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(number);
      record->save_mesh(mesh);
      if(space != NULL)
        record->save_space(space);
      if(sln != NULL)
        record->save_solution(sln);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);

      this->numbered_records.insert(std::pair<unsigned int, CalculationContinuity<Scalar>::Record*>(number, record));
      this->last_record = record;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::add_record(unsigned int number, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> slns, double time_step, double time_step_n_minus_one, double error)
    {
      std::ofstream ofile("onlyNumber.h2d", std::ios_base::app);
      if(ofile)
      {
        ofile << ++this->num << ' ' << number << std::endl;
        ofile.close();
      }
      else
        throw IOCalculationContinuityException(CalculationContinuityException::general, IOCalculationContinuityException::output, "onlyNumber.h2d");

      CalculationContinuity<Scalar>::Record* record = new CalculationContinuity<Scalar>::Record(number);
      record->save_meshes(meshes);
      if(spaces != Hermes::vector<Space<Scalar>*>())
        record->save_spaces(spaces);
      if(slns != Hermes::vector<Solution<Scalar>*>())
        record->save_solutions(slns);
      if(time_step > 0.0)
        record->save_time_step_length(time_step);
      if(time_step_n_minus_one > 0.0)
        record->save_time_step_length_n_minus_one(time_step_n_minus_one);
      if(error > 0.0)
        record->save_error(error);
      this->numbered_records.insert(std::pair<unsigned int, CalculationContinuity<Scalar>::Record*>(number, record));
      this->last_record = record;
    }

    template<typename Scalar>
    CalculationContinuity<Scalar>::Record::Record(double time, unsigned int number) : time(time), number(number)
    {
    }

    template<typename Scalar>
    CalculationContinuity<Scalar>::Record::Record(double time) : time(time), number(0)
    {
    }

    template<typename Scalar>
    CalculationContinuity<Scalar>::Record::Record(unsigned int number) : time(0.0), number(number)
    {
    }

    template<typename Scalar>
    bool CalculationContinuity<Scalar>::have_record_available()
    {
      return this->record_available;
    }

    template<typename Scalar>
    typename CalculationContinuity<Scalar>::Record* CalculationContinuity<Scalar>::get_last_record() const
    {
      if(this->last_record != NULL)
        return this->last_record;
      else
        throw CalculationContinuityException(CalculationContinuityException::general, "no last_record in get_last_record()");
    }

    template<typename Scalar>
    int CalculationContinuity<Scalar>::get_num() const
    {
      return this->num;
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_meshes(Hermes::vector<Mesh*> meshes)
    {
      MeshReaderH2DXML reader;
      for(unsigned int i = 0; i < meshes.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::mesh_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
        try
        {
          reader.save(filename.str().c_str(), meshes[i]);
        }
        catch(std::exception& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::meshes, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
        }
      }
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_mesh(Mesh* mesh)
    {
      MeshReaderH2DXML reader;
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::mesh_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        reader.save(filename.str().c_str(), mesh);
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::meshes, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_spaces(Hermes::vector<Space<Scalar>*> spaces)
    {
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::space_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
        try
        {
          spaces[i]->save(filename.str().c_str());
        }
        catch(std::exception& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
        }
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_space(Space<Scalar>* space)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::space_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        space->save(filename.str().c_str());
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_solutions(Hermes::vector<Solution<Scalar>*> solutions)
    {
      for(unsigned int i = 0; i < solutions.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::solution_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

        try
        {
          solutions[i]->save(filename.str().c_str());
        }
        catch(Hermes::Exceptions::SolutionSaveFailureException& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
        }
      }
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_solution(Solution<Scalar>* solution)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::solution_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        solution->save(filename.str().c_str());
      }
      catch(Hermes::Exceptions::SolutionSaveFailureException& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_time_step_length(double time_step_length_to_save)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::time_step_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ofstream out(filename.str().c_str());
        out << time_step_length_to_save;
        out.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::time_steps, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_time_step_length_n_minus_one(double time_step_length_to_save)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::time_stepNMinusOne_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ofstream out(filename.str().c_str());
        out << time_step_length_to_save;
        out.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::time_steps, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::save_error(double error)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::error_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ofstream out(filename.str().c_str());
        out << error;
        out.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::error, IOCalculationContinuityException::output, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_meshes(Hermes::vector<Mesh*> meshes)
    {
      MeshReaderH2DXML reader;
      for(unsigned int i = 0; i < meshes.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::mesh_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
        try
        {
          reader.load(filename.str().c_str(), meshes[i]);
        }
        catch(Hermes::Exceptions::MeshLoadFailureException& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::meshes, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
      }
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_mesh(Mesh* mesh)
    {
      MeshReaderH2DXML reader;
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::mesh_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

      try
      {
        reader.load(filename.str().c_str(), mesh);
      }
      catch(Hermes::Exceptions::MeshLoadFailureException& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::meshes, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    Hermes::vector<Space<Scalar>*> CalculationContinuity<Scalar>::Record::load_spaces(Hermes::vector<Mesh*> meshes, Hermes::vector<EssentialBCs<Scalar>*> essential_bcs, Hermes::vector<Shapeset*> shapesets)
    {
			Hermes::vector<Space<Scalar>*> spaces;

			if(shapesets == Hermes::vector<Shapeset*>())
        for(unsigned int i = 0; i < meshes.size(); i++)
          shapesets.push_back(NULL);

      for(unsigned int i = 0; i < meshes.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::space_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

        try
        {
          spaces.push_back(Space<Scalar>::load(filename.str().c_str(), meshes[i], false, essential_bcs[i], shapesets[i]));
        }
        catch(Hermes::Exceptions::SpaceLoadFailureException& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
        catch(std::exception& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
      }
    }

    template<typename Scalar>
    Hermes::vector<Space<Scalar>*> CalculationContinuity<Scalar>::Record::load_spaces(Hermes::vector<Mesh*> meshes, Hermes::vector<Shapeset*> shapesets)
    {
			Hermes::vector<Space<Scalar>*> spaces;

      if(shapesets == Hermes::vector<Shapeset*>())
        for(unsigned int i = 0; i < meshes.size(); i++)
          shapesets.push_back(NULL);

      for(unsigned int i = 0; i < meshes.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::space_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

        try
        {
          spaces.push_back(Space<Scalar>::load(filename.str().c_str(), meshes[i], false, NULL, shapesets[i]));
        }
        catch(Hermes::Exceptions::SpaceLoadFailureException& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
        catch(std::exception& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
      }
    }

    template<typename Scalar>
    Space<Scalar>* CalculationContinuity<Scalar>::Record::load_space(Mesh* mesh, EssentialBCs<Scalar>* essential_bcs, Shapeset* shapeset)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::space_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

      try
      {
        return Space<Scalar>::load(filename.str().c_str(), mesh, false, essential_bcs, shapeset);
      }
      catch(Hermes::Exceptions::SpaceLoadFailureException& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::spaces, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_solutions(Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<Space<Scalar>*> spaces)
    {
      if(solutions.size() != spaces.size())
        throw Exceptions::LengthException(1, 2, solutions.size(), spaces.size());
      for(unsigned int i = 0; i < solutions.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::solution_file_name << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
        try
        {
          solutions[i]->load(filename.str().c_str(), spaces[i]);
          solutions[i]->space_type = spaces[i]->get_type();
        }
        catch(Hermes::Exceptions::SolutionLoadFailureException& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
        catch(std::exception& e)
        {
          throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
        }
      }
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_solution(Solution<Scalar>* solution, Space<Scalar>* space)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::solution_file_name << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        solution->load(filename.str().c_str(), space);
        solution->space_type = space->get_type();
      }
      catch(Hermes::Exceptions::SolutionLoadFailureException& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::solutions, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_time_step_length(double & time_step_length)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::time_step_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ifstream in(filename.str().c_str());
        in >> time_step_length;
        in.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::time_steps, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_time_step_length_n_minus_one(double & time_step_length)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::time_stepNMinusOne_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ifstream in(filename.str().c_str());
        in >> time_step_length;
        in.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::time_steps, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    void CalculationContinuity<Scalar>::Record::load_error(double & error)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::error_file_name << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        std::ifstream in(filename.str().c_str());
        in >> error;
        in.close();
      }
      catch(std::exception& e)
      {
        throw IOCalculationContinuityException(CalculationContinuityException::error, IOCalculationContinuityException::input, filename.str().c_str(), e.what());
      }
    }

    template<typename Scalar>
    double CalculationContinuity<Scalar>::Record::get_time()
    {
      return this->time;
    }

    template<typename Scalar>
    unsigned int CalculationContinuity<Scalar>::Record::get_number()
    {
      return this->number;
    }

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::mesh_file_name = "Mesh-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::space_file_name = "Space-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::solution_file_name = "Solution-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::time_step_file_name = "TimeStep_";
    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::time_stepNMinusOne_file_name = "TimeStepNMinusOne_";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::error_file_name = "Error_";

    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_mesh_file_name(std::string mesh_file_nameToSet)
    {
      mesh_file_name = mesh_file_nameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_space_file_name(std::string space_file_nameToSet)
    {
      space_file_name = space_file_nameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_solution_file_name(std::string solution_file_nameToSet)
    {
      solution_file_name = solution_file_nameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_time_step_file_name(std::string time_step_file_nameToSet)
    {
      time_step_file_name = time_step_file_nameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_error_file_name(std::string error_file_nameToSet)
    {
      error_file_name = error_file_nameToSet;
    }

    template class HERMES_API CalculationContinuity<double>;
    template class HERMES_API CalculationContinuity<std::complex<double> >;
  }
}
