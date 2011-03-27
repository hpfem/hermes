#include "hermes2d.h"
#include <GL/glew.h>
#include <GL/freeglut.h>

// This example allows you to try out the automatic zooming capability of ScalarView,
// as well as displaying constant functions and bounding boxes.
// Called as
//    zoom-to-fit <function> <domain>
// it displays a plot of a function specified by the argument <function> over the mesh
// loaded from file <domain>, using several screen sizes. Before going to another screen
// size, the program waits until you finish inspecting the current window and
// close it.
// The argument <function> may attain values from 0 to 3, corresponding to functions that
// are currently implemented in 'test_functions.cpp':
//   0: Constant function:  z = 3.
//   1: Plane parallel to the x-axis with a 50° angle with the xy-plane:  z = tan(50°)*y.
//         After the default rotation of another 40° around the x-axis, a rectangle in the
//         screen plane should be shown (the lighting should be fixed somehow however so that
//         it is filled).
//   2: Cuboid-like function.
//   3: Function with a large range of values.
//         This function is displayed twice for each screen size. Using first the automatic choice of
//         z-axis limits leads to a thin tall plot as ScalarView tries to display all function values
//         (this should be fixed in future by a suitable automatic y-axis scaling according to the value range
//         and/or adjusting the zooming range). The second time, the range is set manually -
//         in future, all above and below the range should be clipped; now, these vertices
//         are still drawn, but using the color corresponding to the set limits (the tester should
//         see it by zooming to the top and bottom faces of the bounding box and noticing how the shading
//         levels vanish into a single shade beyond the box).
//   4: Paraboloid: z = x^2 + y^2.
//         This tests manual setting of limits.
// There are three domains prepared - the square, trapezoid and an l-shape with curved part.
// Functions 2 and 3 don't work with the trapezoid.

#include "definitions.cpp"

extern void init_glut();

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    printf("Please input as this format: zoom-to-fit <function> <domain> \n");
    return -1;
  }

  // Define dimensions of the various tested views.
  init_glut();
  int screen_width = glutGet(GLUT_SCREEN_WIDTH);
  int screen_height = glutGet(GLUT_SCREEN_HEIGHT);
  int rect_size = std::min(screen_height, screen_width);

  int test_dims[][2] = {
    {rect_size / 6, rect_size / 6},
    {rect_size, rect_size},
    {screen_width / 6, screen_height},
    {screen_width / 2, screen_height},
    {screen_width, screen_height / 6},
    {screen_width, screen_height / 2},
  };

  // Window titles for the selected function.
  const std::string titles[5] = {
    "z = c",
    "z = tan(50)*y",
    "Cuboid-like function",
    "Function with a large range",
    "z = x^2 + y^2"
  };
  std::string title;

  bool auto_range = true;      // True to determine the vertical limits from function values.
  double range_min, range_max; // Custom vertical limits.

  // Function for inspection.
  Solution fn;

  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(argv[2], &mesh);

  int fn_id = atoi(argv[1]);
  switch(fn_id) {
    case 0:
      fn.set_exact(&mesh, fn_const);
      title = titles[0];
      break;
    case 1:
      fn.set_exact(&mesh, fn_plane);
      title = titles[1];
      break;
    case 2:
      fn.set_exact(&mesh, fn_cuboid);
      title = titles[2];
      break;
    case 3:
      fn.set_exact(&mesh, fn_bigrange);
      range_min = -2;
      range_max = 4;
      auto_range = false;
      title = titles[3];
      break;
    case 4:
      fn.set_exact(&mesh, fn_paraboloid);
      auto_range = false;
      title = titles[4];
      break;
    default:
      printf("Please set the first argument to a number from 0 to 4: \n");
      for (int i = 0; i < 4; i++)
        printf("%d: %s\n", i, title.c_str());
      return -1;
  }

  if (fn_id == 4) { // Test manual setting bounds for the displayed range.

    ScalarView view(const_cast<char *>(title.c_str()), new WinGeom(screen_width/4, screen_height/4, screen_width/2, screen_width/2));
    view.set_3d_mode(true);

    // Test the behaviour when user enters bigger lower bound.
    view.set_min_max_range(1, 0.5);
    // Show the function.
    view.show(&fn);
    // Wait for the view to be closed.
    View::wait();

    // Test the behaviour when user enters both bounds the same.
    view.set_min_max_range(0.5, 0.5);
    // Show the function.
    view.show(&fn);
    // Wait for the view to be closed.
    View::wait();

  } else { // Test model positioning.

    for (int i = 0; i < 6; i++) {
      ScalarView view(const_cast<char *>(title.c_str()), new WinGeom(0, 0, test_dims[i][0], test_dims[i][1]));
      view.set_3d_mode(true);

      // Show the function.
      view.show(&fn);
      // Wait for the view to be closed.
      View::wait();

      if (!auto_range) {
        char buf[256];
        sprintf(buf, "%s - restricted to (%f,%f)", title.c_str(), range_min, range_max);
        ScalarView view(buf, new WinGeom(0, 0, test_dims[i][0], test_dims[i][1]));
        view.set_min_max_range(range_min, range_max);
        view.set_3d_mode();
        view.show_bounding_box();

        // Show the function.
        view.show(&fn);
        // Wait for the view to be closed.
        View::wait();
      }
    }
  }

  return 0;
}
