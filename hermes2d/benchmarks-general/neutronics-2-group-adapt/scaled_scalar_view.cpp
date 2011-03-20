#include "scaled_scalar_view.h"

void ScaledScalarView::scale(double sc)
{
    /* FIXME: If I uncomment the lines below, I am getting:
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp: In member function ‘void ScaledScalarView::scale(double)’:
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:4: error: ‘class ScaledScalarView’ has no member named ‘lin’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:5: error: ‘class ScaledScalarView’ has no member named ‘mode3d’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:6: error: ‘class ScaledScalarView’ has no member named ‘yscale’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:7: error: ‘class ScaledScalarView’ has no member named ‘contours’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:8: error: ‘class ScaledScalarView’ has no member named ‘cont_step’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:9: error: ‘class ScaledScalarView’ has no member named ‘lin’
/home/certik1/repos/hermes2d/benchmarks/neutronics-2-group-adapt/scaled_scalar_view.cpp:10: error: ‘class ScaledScalarView’ has no member named ‘refresh’


However, the "lin" is declared in src/views/scalar_view.h as a protected member of the class ScalarView, so it should work.
*/

    /*
    this->lin.lock_data();
    if (this->mode3d)
        this->yscale *= sc;
    else if (this->contours)
        this->cont_step *= sc;
    this->lin.unlock_data();
    this->refresh();
    */
}
