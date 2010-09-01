#ifndef __H2D_TEST_CAND_PROJ_H
#define __H2D_TEST_CAND_PROJ_H

/// Test case.
class TestCase {
  int m_func_quad_order; ///< Function order.
  int m_start_quad_order; ///< Start order of an element.
  double** m_poly_matrix; ///< MxN matrix of polygonal coefficients. The value is calculates [y^0, ..., y^(M-1)] * matrix * [x^0, ..., x^(N-1)]^T
public:
  TestCase(int func_quad_order); ///< Contructor. Creates a instance of a test function of a given order.
  ~TestCase(); ///< Destructor. Cleans allocated structures.

  bool should_match(const RefinementSelectors::OptimumSelector::Cand& cand); ///< Returns true if the refinement should match the function.

  std::string title() const; ///< Returns title of the test case. Title is generated automatically.
  int quad_order() const { return m_func_quad_order; } ///< Returns function quad order.
  int start_quad_order() const { return m_start_quad_order; } ///< Returns starting order of a mesh.
  double** poly_matrix() { return m_poly_matrix; }; ///< Returns poly matrix.
};

//test functions
extern scalar func_val(double** poly_matrix, int quad_order, double x, double y); ///< Returns a value f of a polynom f(x,y).
extern scalar func_dx(double** poly_matrix, int quad_order, double x, double y); ///< Returns df(x,y)/dx of a polynom f(x,y).
extern scalar func_dy(double** poly_matrix, int quad_order, double x, double y); ///< Returns df(x,y)/dy of a polynom f(x,y).

#endif
