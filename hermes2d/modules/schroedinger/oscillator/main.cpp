#include "schroedinger.h"

int main(int argc, char* argv[])
{
    ModuleSchroedinger m;

    PotentialHarmonicOscillator p;
    p.set_omega(5);
    m.set_potential(&p);

    Matrix *A, *B;
    //m.assemble(A, B);
}
