#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

#include "exact_solution.cpp"

/* Exact solution */

class CustomExactSolution : public ExactSolutionVector
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionVector(mesh) {};
  ~CustomExactSolution() {};

  scalar2 value(double x, double y) const {
    return  scalar2(x, y);//scalar2(exact0(x, y, dx, dy), exact1(x, y, dx, dy));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {};
  
  virtual Ord ord(Ord x, Ord y) const {
    return Ord(2);
  } 
};
