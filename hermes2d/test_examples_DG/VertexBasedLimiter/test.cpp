#include "definitions.h"

#define SHOW_OUTPUT

// 3.75x^2 + 6.12x + 4.75y^2 + 5.12y + 7.23
class TestPolynomial1 : public ExactSolutionScalar<double>
{
public:
  TestPolynomial1(MeshSharedPtr mesh, int order) : ExactSolutionScalar<double>(mesh), order(order)
  {
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const
  {
    switch(order)
    {
    case 0:
      dx = 0.;
      dy = 0.;
      break;
    case 1:
      dx = 6.12;
      dy = 5.12;
      break;
    case 2:
      dx = 7.5*x + 6.12;
      dy = 9.5*y + 5.12;
      break;
    }
  }

  virtual double value (double x, double y) const
  {
    double result = 7.23;
    if(order > 0)
      result += 6.12 * x + 5.12 * y;
    if(order > 1)
      result += 3.75 * x * x + 4.75 * y * y;

    return result;
  }

  virtual Ord ord(double x, double y) const
  {
    return Ord(5);
  }

  MeshFunction<double>* clone() const
  {
    return new TestPolynomial1(this->mesh, this->order);
  }

  int order;
};

class TestPolynomial2 : public ExactSolutionScalar<double>
{
public:
  TestPolynomial2(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
  {
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const
  {
    if(x < 1)
    {
    if(x > 0.)
    {
      if(y > 0.)
      {
        dx = 2.;
        dy = 2.;
      }
      else
      {
        dx = -2.;
        dy = 2.;
      }
    }
    else
    {
      if(y > 0.)
      {
        dx = 0.;
        dy = 0.;
      }
      else
      {
        dx = 0.;
        dy = 0.;
      }
    }
    }
    else
    {
      if(y > 0.)
      {
        dx = 0.;
        dy = 0.;
      }
      else
      {
        dx = 0.;
        dy = 0.;
      }
    }
  }

  virtual double value (double x, double y) const
  {
    if(x < 1)
    {
    if(x > 0.)
    {
      if(y > 0.)
      {
        return 2.*(x + y);
      }
      else
      {
        return 2.* ((1 - x) + (1 + y));
      }
    }
    else
    {
      if(y > 0.)
      {
        return 3.;
      }
      else
      {
        return 1.;
      }
    }
    }
    else
      if(y > 0.)
      {
        return 5.;
      }
      else
      {
        return -1.;
      }
  }

  virtual Ord ord(double x, double y) const
  {
    return Ord(5);
  }

  MeshFunction<double>* clone() const
  {
    return new TestPolynomial2(this->mesh);
  }
};

class TestPolynomial3 : public ExactSolutionScalar<double>
{
public:
  TestPolynomial3(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh)
  {
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const
  {
    if(x > 1.)
      dx = 3.0;
    else
    {
      if(x > 0.)
        dx = 5.*x;
      else
        dx = 1.;
    }
    dy = 0.;
  }

  virtual double value (double x, double y) const
  {
    if(x > 1.)
      return 3.*x;
    else
    {
      if(x > 0.)
        return 2.5*x*x;
      else
        return x;
    }
  }

  virtual Ord ord(double x, double y) const
  {
    return Ord(5);
  }

  MeshFunction<double>* clone() const
  {
    return new TestPolynomial3(this->mesh);
  }
};

bool test_value(double obtained_value, double expected_value, const char* identifier, double absolute_precision = 1e-6)
{
  if(std::abs(expected_value - obtained_value) < absolute_precision) 
    return true;
  else
  {
    std::cout << "Failed test: " << identifier << std::endl;
    std::cout << "Difference: " << expected_value << " - " << obtained_value << " higher than the precision (" << absolute_precision << ")." << std::endl;
    return false;
  }
}

void test_projections(bool& success)
{
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("test_mesh_projections.xml", mesh);

#ifdef SHOW_OUTPUT
  ScalarView solution_view("Solution", new WinGeom(50, 50, 500, 500));
  Linearizer lin;
#endif

  for(int order = 0; order <= 2; order++)
  {
    SpaceSharedPtr<double> space(new L2Space<double>(mesh, order, new L2ShapesetTaylor));
    double* sln_vector = new double[space->get_num_dofs()];

    MeshFunctionSharedPtr<double>solution(new TestPolynomial1(mesh, order));
    OGProjection<double>::project_global(space, solution, sln_vector);
#ifdef SHOW_OUTPUT
    MeshFunctionSharedPtr<double>projected_solution(new Solution<double>());
    OGProjection<double>::project_global(space, solution, projected_solution);
    solution_view.show(projected_solution);
    lin.save_solution_vtk(projected_solution, "sln.vtk", "sln");
    View::wait_for_keypress();
#endif

    if(order == 0)
    {
      success = test_value(sln_vector[0], 7.23, "projection, order = 0, fn = 0") && success;
    }
    if(order == 1)
    {
      success = test_value(sln_vector[0], 12.85, "projection, order = 1, fn = 0") && success;
      success = test_value(sln_vector[1], 6.12, "projection, order = 1, fn = 1") && success;
      success = test_value(sln_vector[2], 5.12, "projection, order = 1, fn = 2") && success;
    }
    if(order == 2)
    {
      success = test_value(sln_vector[0], 15.6833333333333, "projection, order = 2, fn = 0") && success;
      success = test_value(sln_vector[1], 9.87, "projection, order = 2, fn = 1") && success;
      success = test_value(sln_vector[2], 9.87, "projection, order = 2, fn = 2") && success;
      success = test_value(sln_vector[3], 7.5, "projection, order = 2, fn = 3") && success;
      success = test_value(sln_vector[4], 9.5, "projection, order = 2, fn = 4") && success;
      success = test_value(sln_vector[5], 0.0, "projection, order = 2, fn = 5") && success;
    }
  }
}

void test_centroids_calculation(bool& success, ElementMode2D mode)
{
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  if(mode == HERMES_MODE_TRIANGLE)
    mloader.load("test_mesh_centroids_tris.xml", mesh);
  else
    mloader.load("test_mesh_centroids_quads.xml", mesh);

  for(int order = 0; order <= 2; order++)
  {
    MeshFunctionSharedPtr<double>solution(new TestPolynomial1(mesh, order));
    TestPolynomial1* polynom = dynamic_cast<TestPolynomial1*>(solution.get());

    double x_center, y_center;
    mesh->get_element(0)->get_center(x_center, y_center);
    double val[3];
    val[0] = polynom->value(x_center, y_center);
    polynom->derivatives(x_center, y_center, val[1], val[2]);

    SpaceSharedPtr<double> space(new L2Space<double>(mesh, order, new L2ShapesetTaylor));
    OGProjection<double>::project_global(space, solution, solution);
    Solution<double>* sln = dynamic_cast<Solution<double>*>(solution.get());

    for(int mixed_derivative_index = 0; mixed_derivative_index <=2 ; mixed_derivative_index++)
    {
      if(mode == HERMES_MODE_TRIANGLE)
      {
        double center_value = sln->get_ref_value_transformed(mesh->get_element(0), CENTROID_TRI_X, CENTROID_TRI_Y, 0, mixed_derivative_index);
        success = test_value(center_value, val[mixed_derivative_index], "centroid - triangle") && success;
      }
      else
      {
        double center_value = sln->get_ref_value_transformed(mesh->get_element(0), CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, mixed_derivative_index);
        success = test_value(center_value, val[mixed_derivative_index], "centroid - quad") && success;
      }
    }
  }
}

void test_vertex_extrema_calculation(bool& success)
{
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("test_square.xml", mesh);
  mesh->refine_all_elements();

  MeshFunctionSharedPtr<double>solution_lin(new TestPolynomial2(mesh));
  MeshFunctionSharedPtr<double>solution_quad(new TestPolynomial3(mesh));

#ifdef SHOW_OUTPUT
  ScalarView solution_view_lin("Solution-linear", new WinGeom(50, 50, 500, 500));
  ScalarView solution_view_quad("Solution-quadratic", new WinGeom(600, 50, 500, 500));
  solution_view_lin.show(solution_lin);
  solution_view_quad.show(solution_quad, HERMES_EPS_NORMAL, H2D_FN_DX_0);
  View::wait_for_keypress();
#endif

  SpaceSharedPtr<double> space_lin(new L2Space<double>(mesh, 1, new L2ShapesetTaylor));
  double* sln_vector_lin = new double[space_lin->get_num_dofs()];
  OGProjection<double>::project_global(space_lin, solution_lin, sln_vector_lin);
  MeshFunctionSharedPtr<double>proj_solution_lin(new Solution<double>());
  OGProjection<double>::project_global(space_lin, solution_lin, proj_solution_lin);
  PostProcessing::VertexBasedLimiter limiter_lin(space_lin, sln_vector_lin, 1);
  solution_lin = limiter_lin.get_solution();
  success = test_value(limiter_lin.get_correction_factors().size(), 2, "limiter-quad-size") && success;
  for(int i = 0; i < limiter_lin.get_correction_factors().size(); i++)
  {
#ifdef SHOW_OUTPUT
    std::cout << "Element no.: " << limiter_lin.get_changed_element_ids()[i] << ": " << limiter_lin.get_correction_factors()[i].first << ", " << limiter_lin.get_correction_factors()[i].second << std::endl;
#endif
    success = test_value(limiter_lin.get_correction_factors()[i].second, 0.5, "limiter-lin") && success;
  }

  SpaceSharedPtr<double> space_quad(new L2Space<double>(mesh, 2, new L2ShapesetTaylor));
  double* sln_vector_quad = new double[space_quad->get_num_dofs()];
  OGProjection<double>::project_global(space_quad, solution_quad, sln_vector_quad);
  MeshFunctionSharedPtr<double>proj_solution_quad(new Solution<double>());
  OGProjection<double>::project_global(space_quad, solution_quad, proj_solution_quad);
  PostProcessing::VertexBasedLimiter limiter_quad(space_quad, sln_vector_quad, 2);
#ifdef SHOW_OUTPUT
  limiter_quad.set_verbose_output(true);
#endif
  solution_quad = limiter_quad.get_solution();
  success = test_value(limiter_quad.get_correction_factors().size(), 2, "limiter-quad-size") && success;
  for(int i = 0; i < limiter_quad.get_correction_factors().size(); i++)
  {
#ifdef SHOW_OUTPUT
    std::cout << "Element no.: " << limiter_quad.get_changed_element_ids()[i] << ": " << limiter_quad.get_correction_factors()[i].first << ", " << limiter_quad.get_correction_factors()[i].second << std::endl;
#endif
    success = test_value(limiter_quad.get_correction_factors()[i].second, 0.2, "limiter-quad") && success;
  }

#ifdef SHOW_OUTPUT
  solution_view_lin.show(proj_solution_lin);
  solution_view_quad.show(proj_solution_quad, HERMES_EPS_NORMAL, H2D_FN_DX_0);
  View::wait_for_keypress();
  
  solution_view_lin.show(solution_lin);
  solution_view_quad.show(solution_quad, HERMES_EPS_NORMAL, H2D_FN_DX_0);
  View::wait_for_keypress();
#endif
}

int test()
{
  bool success = true;
  test_projections(success);
  test_centroids_calculation(success, HERMES_MODE_TRIANGLE);
  test_centroids_calculation(success, HERMES_MODE_QUAD);
  test_vertex_extrema_calculation(success);

  if(success)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }
}