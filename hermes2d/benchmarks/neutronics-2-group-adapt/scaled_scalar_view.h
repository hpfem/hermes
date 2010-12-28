#include "hermes2d.h"

class HERMES_API ScaledScalarView : public ScalarView
{
public:
    ScaledScalarView(char* title = "ScaledScalarView", WinGeom *wg = NULL)
                  : ScalarView(title, wg) {
    }

    void scale(double sc = 1e3);
};
