#include <hermes2d.h>

using namespace RefinementSelectors;

/** \addtogroup t_bubbles Bubbles
 *  \{
 *  \brief This test tests bubble functions on quads and their orders.
 *
 *  The test succeeds if:
 *   - There are defined bubbles for every permutation of orders from a given range.
 *   - (Integration) orders of retrieved bubbles does not fit in a given range of orders.
 *
 *  Output table legend:
 *   - '-': Test was not done.
 *   - ' ': Test succesfull.
 *   - 'D': Bubbles are not defined for a given range of orders
 *   - 'S': (Integration) orders of retrieved bubbles does not fit in a given range of orders.
 */

/* global definitions */
#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS 0 ///< Test return code if success.
#define ERROR_FAILURE -1 ///< Test return code if fails.
#define delete_not_null(__ptr) if (__ptr != NULL) delete __ptr; ///< Deletes an instance if the pointer is not NULL.

#define H2D_TEST_NOT_DONE -1 ///< Flag: Combination of orders was not tested.
#define H2D_TEST_SUCCESS 0x00 ///< Flag: Combination of orders was tested successfully.
#define H2D_NOT_DEFINED 0x01  ///< Flag: A given combination is not defined
#define H2D_INVALID_SHAPE 0x02  ///< Flag: Invalid shape function.

/* global variables */
std::string space_name; ///< A name of the space.
Shapeset* shapeset = NULL; ///< A shapeset used by the test.
Range<int> order_range; ///< A order range in which bubbles are examined
int intr_order_shift; ///< A shift of integration order from the order of the element.

/// Cleanup
void cleanup() {
  delete_not_null(shapeset);
}

/// Initialize internal data structures: H1 space.
bool init_h1() {
  shapeset = new H1Shapeset();
  intr_order_shift = 0;
  order_range = Range<int>(2, shapeset->get_max_order());
  space_name = "H1";
  return true;
}

/// Initialize internal data structures: L2 space.
bool init_l2() {
  shapeset = new L2Shapeset();
  intr_order_shift = 0;
  order_range = Range<int>(2, shapeset->get_max_order());
  space_name = "L2";
  return true;
}

/// Initialize internal data structures: Hcurl space.
bool init_hcurl() {
  shapeset = new HcurlShapeset();
  intr_order_shift = 1;
  order_range = Range<int>(1, shapeset->get_max_order() - intr_order_shift);
  space_name = "Hcurl";
  return true;
}

/// Prints failure matrix
void show_fail_matrix(int** fail_matrix, const std::string& space_name, const int max_quad_order) {
  const int max_order_h = H2D_GET_H_ORDER(max_quad_order), max_order_v = H2D_GET_V_ORDER(max_quad_order);
#define NUMBER_W 2 ///< A size of a cell in the table.
#define NUMBER_FMT "%02d" ///< A format for writing orders.
  char buff_number[1024];

  info("!Test summary (V/H): %s", space_name.c_str());

  //print header
  {
    std::stringstream str;
    for(int i = 0; i < NUMBER_W; i++)
      str << ' ';
    for(int i = 0; i <= max_order_h; i++) {
      str << '|';
      sprintf(buff_number, NUMBER_FMT, i);
      str << buff_number;
    }
    info(" %s", str.str().c_str());
  }

  //print body
  for(int i = 0; i <= max_order_v; i++) {
    //build row head
    std::stringstream str;
    sprintf(buff_number, NUMBER_FMT, i);
    str << buff_number;

    //build row body
    for(int k = 0; k <= max_order_h; k++) {
      str << '|';
      if (fail_matrix[i][k] == H2D_TEST_NOT_DONE) {
        for(int j = 0; j < NUMBER_W; j++)
          str << '-';
      }
      else {
        if ((fail_matrix[i][k] & H2D_NOT_DEFINED) != 0)
          str << 'D';
        else
          str << ' ';
        if ((fail_matrix[i][k] & H2D_INVALID_SHAPE) != 0)
          str << 'S';
        else
          str << ' ';
        for(int j = 2; j < NUMBER_W; j++)
          str << ' ';
      }
    }

    //print row
    info(" %s", str.str().c_str());
  }
}

/// Test
bool test() {
  bool failed = false;
  shapeset->set_mode(H2D_MODE_QUAD);

  //report name
  info("!Space: %s", space_name.c_str());

  //prepare cases
  int min_quad_order = -1, max_quad_order = -1;
  min_quad_order = H2D_MAKE_QUAD_ORDER(order_range.lower(), order_range.lower());
  max_quad_order = H2D_MAKE_QUAD_ORDER(order_range.upper(), order_range.upper());
  OrderPermutator order_perm(min_quad_order, max_quad_order, false);

  //prepare place for result summary
  int end_order_v = H2D_GET_V_ORDER(order_perm.get_end_quad_order());
  int end_order_h = H2D_GET_H_ORDER(order_perm.get_end_quad_order());
  int** fail_matrix = new_matrix<int>(end_order_v+1, end_order_h+1);
  for(int i = 0; i <= end_order_v; i++)
    for(int k = 0; k <= end_order_h; k++)
      fail_matrix[i][k] = H2D_TEST_NOT_DONE;

  //process cases
  do {
    int test_result = H2D_TEST_SUCCESS;

    //retrieve all bubbles
    const int num_indices = shapeset->get_num_bubbles(order_perm.get_quad_order());
    const int* indices = shapeset->get_bubble_indices(order_perm.get_quad_order());
    if (num_indices == 0 || indices == NULL)
      test_result = H2D_NOT_DEFINED;
    else { //bubbles are defined: inspect them
      for(int i = 0; i < num_indices; i++) {
        //retrieve order of a shapeset
        int shape_inx = indices[i];
        int shape_quad_order = shapeset->get_order(shape_inx);
        int shape_order_h = H2D_GET_H_ORDER(shape_quad_order), shape_order_v = H2D_GET_V_ORDER(shape_quad_order);

        //check if shapeset is valid
        if (shape_order_h > order_perm.get_order_h() + intr_order_shift || shape_order_v > order_perm.get_order_v() + intr_order_shift) {
          //report
          verbose("Bubbles of an element order (%d/%d) requested but bubble #%d of order (%d/%d) retrieved.", order_perm.get_order_h(), order_perm.get_order_v(), shape_inx, shape_order_h, shape_order_v);

          //set test result
          test_result = H2D_INVALID_SHAPE;
        }
      }
    }

    if (test_result != H2D_TEST_SUCCESS)
      failed = true;

    fail_matrix[order_perm.get_order_v()][order_perm.get_order_h()] = test_result;
  } while(order_perm.next());

  //print result matrix
  show_fail_matrix(fail_matrix, space_name, order_perm.get_end_quad_order());

  //clenup
  delete[] fail_matrix;

  return !failed;
}

/// Test entry-point.
int main(int argc, char* argv[]) {
  bool test_success = true;

  //H1
  if (init_h1())
    test_success &= test();
  cleanup();
  if (!test_success)
    goto quit;

  //L2
  if (init_l2())
    test_success &= test();
  cleanup();
  if (!test_success)
    goto quit;

  //Hcurl
  if (init_hcurl())
    test_success &= test();
  cleanup();
  if (!test_success)
    goto quit;

quit:
  if (test_success)
  {
    info("!Test: Success");
    return ERROR_SUCCESS;
  }
  else {
    info("!Test: Failed!");
    return ERROR_FAILURE;
  }
}

/// \}
