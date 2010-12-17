#include "hermes2d.h"
#include "schroedinger.h"

int main(int argc, char* argv[])
{
    ModuleSchroedinger m;

    RCP<PotentialHarmonicOscillator> p = rcp(new PotentialHarmonicOscillator());
    p->set_omega(5);
    m.set_potential(p);

    RCP<CSCMatrix> A = rcp(new CSCMatrix());
    RCP<CSCMatrix> B = rcp(new CSCMatrix());
    m.assemble(A, B);
    EigenSolver es(A, B);
    es.solve();
    es.print_eigenvalues();
}
