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
            ifile >> last_time >> last_number;
            record_available = true;
            break;
          case onlyTime:
            ifile >> last_time;
            record_available = true;
            break;
          case onlyNumber:
            ifile >> last_number;
            record_available = true;
            break;
          }
          num++;
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
        ofile << time << ' ' << number << std::endl;
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
        ofile << time << ' ' << number << std::endl;
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
        ofile << time << std::endl;
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
        ofile << time << std::endl;
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
        ofile << number << std::endl;
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
        ofile << number << std::endl;
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
        filename << CalculationContinuity<Scalar>::meshFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::meshFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
        filename << CalculationContinuity<Scalar>::spaceFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::spaceFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
        filename << CalculationContinuity<Scalar>::solutionFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

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
      filename << CalculationContinuity<Scalar>::solutionFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::timeStepFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::timeStepNMinusOneFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::errorFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
        filename << CalculationContinuity<Scalar>::meshFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::meshFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

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
    void CalculationContinuity<Scalar>::Record::load_spaces(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<SpaceType> space_types, Hermes::vector<Mesh*> meshes, Hermes::vector<EssentialBCs<Scalar>*> essential_bcs, Hermes::vector<Shapeset*> shapesets)
    {
      if(shapesets == Hermes::vector<Shapeset*>())
        for(unsigned int i = 0; i < spaces.size(); i++)
          shapesets.push_back(NULL);

      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::spaceFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

        spaces[i]->free();
        try
        {
          switch(space_types[i])
          {
          case HERMES_H1_SPACE:
            dynamic_cast<H1Space<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], essential_bcs[i], shapesets[i]);
            break;
          case HERMES_HCURL_SPACE:
            dynamic_cast<HcurlSpace<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], essential_bcs[i], shapesets[i]);
            break;
          case HERMES_HDIV_SPACE:
            dynamic_cast<HdivSpace<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], essential_bcs[i], shapesets[i]);
            break;
          case HERMES_L2_SPACE:
            dynamic_cast<L2Space<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], shapesets[i]);
            break;
          }
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
    void CalculationContinuity<Scalar>::Record::load_spaces(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<SpaceType> space_types, Hermes::vector<Mesh*> meshes, Hermes::vector<Shapeset*> shapesets)
    {
      if(shapesets == Hermes::vector<Shapeset*>())
        for(unsigned int i = 0; i < spaces.size(); i++)
          shapesets.push_back(NULL);

      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        std::stringstream filename;
        filename << CalculationContinuity<Scalar>::spaceFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

        spaces[i]->free();

        try
        {
          switch(space_types[i])
          {
          case HERMES_H1_SPACE:
            dynamic_cast<H1Space<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], shapesets[i]);
            break;
          case HERMES_HCURL_SPACE:
            dynamic_cast<HcurlSpace<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], shapesets[i]);
            break;
          case HERMES_HDIV_SPACE:
            dynamic_cast<HdivSpace<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], shapesets[i]);
            break;
          case HERMES_L2_SPACE:
            dynamic_cast<L2Space<Scalar>*>(spaces[i])->load(filename.str().c_str(), meshes[i], shapesets[i]);
            break;
          }
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
    void CalculationContinuity<Scalar>::Record::load_space(Space<Scalar>* space, SpaceType space_type, Mesh* mesh, EssentialBCs<Scalar>* essential_bcs, Shapeset* shapeset)
    {
      std::stringstream filename;
      filename << CalculationContinuity<Scalar>::spaceFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";

      space->free();

      try
      {
        switch(space_type)
        {
        case HERMES_H1_SPACE:
          if(essential_bcs == NULL)
            dynamic_cast<H1Space<Scalar>*>(space)->load(filename.str().c_str(), mesh, shapeset);
          else
            dynamic_cast<H1Space<Scalar>*>(space)->load(filename.str().c_str(), mesh, essential_bcs, shapeset);
          break;
        case HERMES_HCURL_SPACE:
          if(essential_bcs == NULL)
            dynamic_cast<HcurlSpace<Scalar>*>(space)->load(filename.str().c_str(), mesh, shapeset);
          else
            dynamic_cast<HcurlSpace<Scalar>*>(space)->load(filename.str().c_str(), mesh, essential_bcs, shapeset);
          break;
        case HERMES_HDIV_SPACE:
          if(essential_bcs == NULL)
            dynamic_cast<HdivSpace<Scalar>*>(space)->load(filename.str().c_str(), mesh, shapeset);
          else
            dynamic_cast<HdivSpace<Scalar>*>(space)->load(filename.str().c_str(), mesh, essential_bcs, shapeset);
          break;
        case HERMES_L2_SPACE:
          dynamic_cast<L2Space<Scalar>*>(space)->load(filename.str().c_str(), mesh, shapeset);
          break;
        }
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
        filename << CalculationContinuity<Scalar>::solutionFileName << i << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
        try
        {
          solutions[i]->load(filename.str().c_str(), spaces[i]->get_mesh());
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
      filename << CalculationContinuity<Scalar>::solutionFileName << 0 << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
      try
      {
        solution->load(filename.str().c_str(), space->get_mesh());
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
      filename << CalculationContinuity<Scalar>::timeStepFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::timeStepNMinusOneFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
      filename << CalculationContinuity<Scalar>::errorFileName << '_' << (std::string)"t = " << this->time << (std::string)"n = " << this->number << (std::string)".h2d";
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
    std::string CalculationContinuity<Scalar>::meshFileName = "Mesh-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::spaceFileName = "Space-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::solutionFileName = "Solution-";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::timeStepFileName = "TimeStep_";
    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::timeStepNMinusOneFileName = "TimeStepNMinusOne_";

    template<typename Scalar>
    std::string CalculationContinuity<Scalar>::errorFileName = "Error_";

    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_meshFileName(std::string meshFileNameToSet)
    {
      meshFileName = meshFileNameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_spaceFileName(std::string spaceFileNameToSet)
    {
      spaceFileName = spaceFileNameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_solutionFileName(std::string solutionFileNameToSet)
    {
      solutionFileName = solutionFileNameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_timeStepFileName(std::string timeStepFileNameToSet)
    {
      timeStepFileName = timeStepFileNameToSet;
    }
    template<typename Scalar>
    void CalculationContinuity<Scalar>::set_errorFileName(std::string errorFileNameToSet)
    {
      errorFileName = errorFileNameToSet;
    }

    template class HERMES_API CalculationContinuity<double>;
    template class HERMES_API CalculationContinuity<std::complex<double> >;
  }
}