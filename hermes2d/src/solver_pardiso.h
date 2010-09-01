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

#ifndef __H2D_SOLVER_PARDISO_H
#define __H2D_SOLVER_PARDISO_H

#ifdef AIX
#define F77_FUNC(func)  func
#else
#define F77_FUNC(func)  func ## _
#endif


// PARDISO prototype
extern  int F77_FUNC(pardisoinit)
    (void *, int *, int *);

extern  int F77_FUNC(pardiso)
    (void *, int *, int *, int *, int *, int *,
     double *, int *, int *, int *, int *, int *,
     int *, double *, double *, int *);


//// PardisoSolver /////////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "solver.h"


/// \brief Pardiso solver wrapper class
///
class PardisoSolver : public Solver
{
public:

  PardisoSolver() { msglvl = 0;
                    memset(default_iparm, 0, sizeof(default_iparm)); }

  void set_iparm(int idx, double val) { default_iparm[idx] = val; }
  void set_msglvl(int lvl) { msglvl = lvl; }


protected:

  virtual bool is_row_oriented()  { return true; }
  virtual bool handles_symmetry() { return true; }

  int msglvl;
  int default_iparm[64]; // TODO

  struct Data
  {
    void* pt[64];
    int iparm[64];
    int mtype;
    int* Ap1;
    int* Ai1;
    int ndofs;
  };


  virtual void* new_context(bool sym)
  {
    Data* data = new Data;
    memset(data->pt, 0, sizeof(data->pt));
    data->Ap1 = NULL;
    data->Ai1 = NULL;

    // matrix type (structurally symmetric / unsymmetric)
    #ifndef H2D_COMPLEX
    data->mtype = sym ? 1 : 11;
    #else
    data->mtype = sym ? 3 : 13;
    #endif

    F77_FUNC(pardisoinit)(data->pt, &data->mtype, data->iparm);

    // numbers of processors, value of OMP_NUM_THREADS
    char* var = getenv("OMP_NUM_THREADS");
    data->iparm[2] = var ? atoi(var) : 1;

    return data;
  }

  virtual void free_context(void* ctx)
  {
    delete (Data*) ctx;
  }


  virtual bool analyze(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym)
  {
    Data* data = (Data*) ctx;
    if (data->Ap1 != NULL) free_data(ctx);

    // convert Ap and Ai to Fortran one-based indexing - this sucks but there's no other way
    data->Ap1 = new int[n+1];
    data->Ai1 = new int[Ap[n]];
    if (data->Ap1 == NULL || data->Ai1 == NULL) error("Out of memory.");
    for (int i = 0; i <= n; i++)
      data->Ap1[i] = Ap[i] + 1;
    for (int i = 0; i < Ap[n]; i++)
      data->Ai1[i] = Ai[i] + 1;
    data->ndofs = n;

    // perform symbolic analysis
    int phase = 11, maxfct = 1, mnum = 1, idum = 0, nrhs = 1, error;
    double ddum = 0.0;
    verbose("Pardiso: analyzing matrix...");
    F77_FUNC(pardiso)(data->pt, &maxfct, &mnum, &data->mtype, &phase,
                      &n, Ax, data->Ap1, data->Ai1,
                      &idum, &nrhs, data->iparm, &msglvl, &ddum, &ddum, &error);
    print_status(error);
    return error == 0;
  }


  virtual bool factorize(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym)
  {
    Data* data = (Data*) ctx;
    int phase = 22, maxfct = 1, mnum = 1, idum = 0, nrhs = 1, error;
    double ddum = 0.0;
    verbose("Pardiso: factorizing matrix...");
    F77_FUNC(pardiso)(data->pt, &maxfct, &mnum, &data->mtype, &phase,
                      &n, Ax, data->Ap1, data->Ai1,
                      &idum, &nrhs, data->iparm, &msglvl, &ddum, &ddum, &error);
    print_status(error);
    return error == 0;
  }


  virtual bool solve(void* ctx, int n, int* Ap, int* Ai, scalar* Ax, bool sym,
                     scalar* RHS, scalar* vec)
  {
    Data* data = (Data*) ctx;
    int phase = 33, maxfct = 1, mnum = 1, idum = 0, nrhs = 1, error;
    double ddum = 0.0;
    verbose("Pardiso: solving system...");
    F77_FUNC(pardiso)(data->pt, &maxfct, &mnum, &data->mtype, &phase,
                      &n, Ax, data->Ap1, data->Ai1,
                      &idum, &nrhs, data->iparm, &msglvl, RHS, vec, &error);
    print_status(error);
    return error == 0;
  }


  virtual void free_data(void* ctx)
  {
    Data* data = (Data*) ctx;
    if (data->Ap1 || data->Ai1)
    {
      // release factorization data
      int phase = -1, maxfct = 1, mnum = 1, idum = 0, nrhs = 1, error;
      double ddum = 0.0;
      F77_FUNC(pardiso)(data->pt, &maxfct, &mnum, &data->mtype, &phase,
                        &data->ndofs, &ddum, data->Ap1, data->Ai1,
                        &idum, &nrhs, data->iparm, &msglvl, &ddum, &ddum, &error);

      // release one-based arrays
      if (data->Ap1 != NULL) { delete [] data->Ap1;  data->Ap1 = NULL; }
      if (data->Ai1 != NULL) { delete [] data->Ai1;  data->Ai1 = NULL; }
    }
  }


  virtual void print_status(int error)
  {
    switch (error)
    {
      case 0:  return;
      case -1: error("Input inconsistent.");
      case -2: error("Not enough memory.");
      case -3: error("Reordering problem.");
      case -4: error("Zero pivot, numerical fact. or iterative refinement problem.");
      case -5: error("Internal error.");
      case -6: error("Preordering failed.");
      case -7: error("Diagonal matrix problem.");
      case -8: error("32-bit overflow problem.");
      default: error("Uknown error.");
    }
  }


};



#undef F77_FUNC

#endif
