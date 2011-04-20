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

#include "../h2d_common.h"
#include "../solution.h"
#include "../discrete_problem.h"
#include "../refmap.h"
#include "../shapeset_h1_all.h"
#include "../quad_all.h"
#include "../h1.h"
#include "../matrix_old.h"
#include "adapt_h1.h"
#include "adapt_ortho_h1.h"
#include "../traverse.h"
#include "../norm.h"
#include "../element_to_refine.h"
#include "../ref_selectors/selector.h"
#include "../ref_selectors/h1_uniform_hp.h"

using namespace std;

H1OrthoHP::H1OrthoHP(int num, Space* space1, Space* space2, Space* space3, Space* space4, Space* space5,
                     Space* space6, Space* space7, Space* space8, Space* space9, Space* space10)
                     : H1AdaptHP(num, space1, space2, space3, space4, space5, space6, space7, space8, space9, space10) {};


bool H1OrthoHP::adapt(double thr, int strat, int adapt_type, bool iso_only, int regularize,
                      double conv_exp, int max_order, bool same_orders, double to_be_processed)
{
  //create refinement oracle
  RefinementSelectors::Selector* refinement_selector = NULL;
  if (adapt_type == 2)
    refinement_selector = new RefinementSelectors::H1OnlyP(max_order);
  else if (adapt_type == 1 && iso_only)
    refinement_selector = new RefinementSelectors::H1OnlyH();
  else
    refinement_selector = new RefinementSelectors::H1UniformHP(iso_only, RefinementSelectors::H2DRS_CAND_HP, conv_exp, max_order);

  //call adaptivity
  bool result = H1AdaptHP::adapt(thr, strat, refinement_selector, regularize, same_orders, to_be_processed);

  //cleanup
  delete refinement_selector;

  return result;
}
