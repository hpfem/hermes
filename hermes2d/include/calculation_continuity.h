// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file calculation_CalculationContinuity.h
\brief Calculation CalculationContinuity functionality.
*/

#include "config.h"
#include "compat.h"
#include "function/solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Own exception class to catch all potential exceptions that might occur in saving / loading entities using the subsequent classes.
    class HERMES_API CalculationContinuityException : public Hermes::Exceptions::Exception
    {
    public:
      enum exceptionEntityType
      {
        meshes,
        spaces,
        solutions,
        time_steps,
        error,
        general
      };
      CalculationContinuityException();
      CalculationContinuityException(exceptionEntityType type, const char * msg);
      void init(exceptionEntityType type, const char * msg);
    };

    /// A derived exception for I/O
    class HERMES_API IOCalculationContinuityException : public CalculationContinuityException
    {
    public:
      enum inputOutput
      {
        input,
        output
      };
      IOCalculationContinuityException(exceptionEntityType type, inputOutput inputOutput, const char * filename);
      IOCalculationContinuityException(exceptionEntityType type, inputOutput inputOutput, const char * filename, const char * reason);
    };

    /// Class used for resuming an interrupted calculation.
    /// Its purpose is to store everything necessary to resume it from a certain point.
    template<typename Scalar>
    class HERMES_API CalculationContinuity
    {
    public:
      /// Choose an identification method of records.
      /// Either both per time step and per a number, or just by one of these.
      enum IdentificationMethod
      {
        timeAndNumber,
        onlyTime,
        onlyNumber
      };

      CalculationContinuity(IdentificationMethod identification_method);

      /// One record of the calculation. Stores every information to resume a calculation from this one point.
      class HERMES_API Record
      {
      public:
        /// Constructors.
        Record(double time, unsigned int number);
        Record(double time);
        Record(unsigned int number);

        /// Saves vector of meshes.
        void save_meshes(Hermes::vector<Mesh*> meshes);
        /// Saves one mesh.
        void save_mesh(Mesh* mesh);

        /// Saves vector of spaces.
        void save_spaces(Hermes::vector<Space<Scalar>*> spaces);
        /// Saves one space.
        void save_space(Space<Scalar>* space);

        /// Saves vector of solutions.
        void save_solutions(Hermes::vector<Solution<Scalar>*> solutions);
        /// Saves one solution.
        void save_solution(Solution<Scalar>* solution);

        /// Saves the time step length.
        void save_time_step_length(double time_step_length_to_save);
        void save_time_step_length_n_minus_one(double time_step_length_to_save);

        /// Saves the spatial error estimate.
        void save_error(double error);

        /// Loads vector of meshes.
        void load_meshes(Hermes::vector<Mesh*> meshes);
        /// Loads one mesh.
        void load_mesh(Mesh* mesh);

        /// Loads vector of spaces.
        void load_spaces(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<SpaceType> space_types, Hermes::vector<Mesh*> meshes, Hermes::vector<EssentialBCs<Scalar>*> essential_bcs, Hermes::vector<Shapeset*> shapeset = Hermes::vector<Shapeset*>());

        /// Loads vector of spaces.
        /// Version without essential BCs.
        void load_spaces(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<SpaceType> space_types, Hermes::vector<Mesh*> meshes, Hermes::vector<Shapeset*> shapeset = Hermes::vector<Shapeset*>());

        /// Loads one space.
        void load_space(Space<Scalar>* space, SpaceType space_type, Mesh* mesh, EssentialBCs<Scalar>* essential_bcs = NULL, Shapeset* shapeset = NULL);

        /// Loads vector of solutions.
        void load_solutions(Hermes::vector<Solution<Scalar>*> solutions, Hermes::vector<Space<Scalar>*> spaces);
        /// Loads one solution.
        void load_solution(Solution<Scalar>* solution, Space<Scalar>* space);

        /// Loads the time step length.
        void load_time_step_length(double & time_step_length);
        void load_time_step_length_n_minus_one(double & time_step_length);

        /// Loads the spatial error estimate.
        void load_error(double & error);

        /// Returns time.
        double get_time();

        /// Returns the number of the current record.
        unsigned int get_number();

      private:
        /// Storage of filenames of needed mesh files.
        Hermes::vector<std::string> meshFiles;
        /// Storage of filenames of needed space files.
        Hermes::vector<std::string> spaceFiles;
        /// Storage of filenames of needed solution files.
        Hermes::vector<std::string> solutionFiles;

        /// Internals. Used for identifying.
        double time;
        unsigned int number;
      };

      /// Add a record.
      /// See records.
      void add_record(double time, unsigned int number, Mesh* mesh, Space<Scalar>* space = NULL, Solution<Scalar>* sln = NULL, double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);
      void add_record(double time, unsigned int number, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces = Hermes::vector<Space<Scalar>*>(), Hermes::vector<Solution<Scalar>*> slns = Hermes::vector<Solution<Scalar>*>(), double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);

      /// Add a record.
      /// See time_records.
      void add_record(double time, Mesh* mesh, Space<Scalar>* space = NULL, Solution<Scalar>* sln = NULL, double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);
      void add_record(double time, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces = Hermes::vector<Space<Scalar>*>(), Hermes::vector<Solution<Scalar>*> slns = Hermes::vector<Solution<Scalar>*>(), double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);

      /// Add a record.
      /// See numbered_records.
      void add_record(unsigned int number, Mesh* mesh, Space<Scalar>* space = NULL, Solution<Scalar>* sln = NULL, double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);
      void add_record(unsigned int number, Hermes::vector<Mesh*> meshes, Hermes::vector<Space<Scalar>*> spaces = Hermes::vector<Space<Scalar>*>(), Hermes::vector<Solution<Scalar>*> slns = Hermes::vector<Solution<Scalar>*>(), double time_step = 0.0, double time_step_n_minus_one = 0.0, double error = 0.0);

      /// Returns the value of record_available.
      /// See record_available.
      bool have_record_available();

      /// Returns a pointer to the last record.
      Record* get_last_record() const;

      /// Returns the count of records.
      int get_num() const;

    private:
      /// Names for the file stored.
      static std::string meshFileName;
      static std::string spaceFileName;
      static std::string solutionFileName;
      static std::string timeStepFileName;
      static std::string timeStepNMinusOneFileName;
      static std::string errorFileName;

      /// Setting of the names for the file stored.
      static void set_meshFileName(std::string meshFileNameToSet);
      static void set_spaceFileName(std::string spaceFileNameToSet);
      static void set_solutionFileName(std::string solutionFileNameToSet);
      static void set_timeStepFileName(std::string timeStepFileNameToSet);
      static void set_errorFileName(std::string errorFileNameToSet);

      /// For time dependent adaptive problems.
      std::map<std::pair<double, unsigned int>, Record*> records;

      /// Just for time dependent problems.
      std::map<double, Record*> time_records;

      /// Possibly for adaptive solution of elliptic problems.
      std::map<unsigned int, Record*> numbered_records;

      /// Last added record.
      Record* last_record;

      /// If there is an old record on the disk.
      bool record_available;

      /// Determines the identification method.
      IdentificationMethod identification_method;

      /// Count of records.
      int num;

      friend class Record;
    };
  }
}