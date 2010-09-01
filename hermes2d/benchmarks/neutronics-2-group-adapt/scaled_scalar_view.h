#include "hermes2d.h"

class H2D_API ScaledScalarView : public ScalarView
{
public:
    ScaledScalarView(const char* title = "ScaledScalarView",
            DEFAULT_WINDOW_POS): ScalarView(title, x, y, width, height) {
    }

    void scale(double sc = 1e3);
};
