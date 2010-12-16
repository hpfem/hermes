#include "schroedinger.h"
#include "matrix_csc.h"

int main(int argc, char* argv[])
{
    ModuleSchroedinger m;

    PotentialHarmonicOscillator p;
    p.set_omega(5);
    m.set_potential(&p);

    CSCMatrix *A, *B;
    //m.assemble(A, B);
}
