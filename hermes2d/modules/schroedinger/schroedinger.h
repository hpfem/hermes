#ifndef __H2D_SCHROEDINGER_H
#define __H2D_SCHROEDINGER_H

#include "hermes2d.h"

class Potential {
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

class PotentialHarmonicOscillator: public Potential {
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

class HERMES_API ModuleSchroedinger {
public:
    ModuleSchroedinger() {
        this->potential = NULL;
    }

    ~ModuleSchroedinger() {
    }

    void set_potential(Potential *potential) {
        this->potential = potential;
    }

    // Solve the problem and return the solution.
    void assemble(const Matrix &A, const Matrix &B);

private:
    Potential *potential;
};

#endif
