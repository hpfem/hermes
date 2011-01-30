#ifndef __H2D_SCHROEDINGER_H
#define __H2D_SCHROEDINGER_H

// The module is a library.
#ifndef EXPORT_HERMES_MODULE
#define EXPORT_HERMES_MODULE
#endif

#include "config.h"

#include "hermes2d.h"

namespace Schroedinger {

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;

class HERMES_MODULE_API Potential {
    public:
        // Returns the value in the point "x" and "y"
        virtual double get_value(double x, double y) = 0;

        // Returns the values in the arrays "x" and "y". You can override this
        // method if you can make this faster for your particular problem.
        // Otherwise it just calls get_value() in a loop.
        virtual void get_values(int n, double *x, double *y, double *values) {
            for (int i=0; i < n; i++)
                values[i] = this->get_value(x[i], y[i]);
        }
};

class HERMES_MODULE_API PotentialHarmonicOscillator: public Potential {
    public:
        PotentialHarmonicOscillator() {
            this->omega = 1.0;
        }

        double get_value(double x, double y) {
            return 0.5 * this->omega * (x*x + y*y);
        }

        void set_omega(double omega) {
            this->omega = omega;
        }
    private:
        double omega;
};

class HERMES_MODULE_API ModuleSchroedinger {
public:
    ModuleSchroedinger() {
        this->potential = null;
    }

    ~ModuleSchroedinger() {
    }

    void set_potential(const RCP<Potential> &potential) {
        this->potential = potential;
    }

    void assemble(const Ptr<SparseMatrix> &A, const Ptr<SparseMatrix> &B);

private:
    RCP<Potential> potential;
};

}

#endif
