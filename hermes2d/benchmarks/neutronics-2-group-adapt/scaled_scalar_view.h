#include "hermes2d.h"

class HERMES_API ScaledScalarView : public ScalarView
{
public:
    ScaledScalarView(const char* title = "ScaledScalarView",
		     DEFAULT_WINDOW_POS): ScalarView(title, new WinGeom(x, y, width, height)) {
    }

    void scale(double sc = 1e3);
};
