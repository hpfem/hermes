#include "weakform/weakform.h"
#include <integrals/integrals_h1.h>
#include "boundaryconditions/essential_bcs.h"

class EssentialBCNonConst : public EssentialBC {
public:
    EssentialBCNonConst(std::string marker) : EssentialBC(Hermes::vector<std::string>())
    {
        markers.push_back(marker);
    }

    ~EssentialBCNonConst() {};

    inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

    virtual scalar value(double x, double y) const
    {
        return 100*cos(y*M_PI/0.1);
    }
};

class WeakFormHelmholtz : public WeakForm
{
public:
    // Problems parameters
    WeakFormHelmholtz(double eps, double mu, double omega, double sigma, double beta, double E0, double h) : WeakForm(2)
    {
        add_matrix_form(new MatrixFormHelmholtzEquation_real_real(0, 0, eps, mu));
        add_matrix_form(new MatrixFormHelmholtzEquation_real_imag(0, 1, mu, omega, sigma));
        add_matrix_form(new MatrixFormHelmholtzEquation_imag_real(1, 0, mu, omega, sigma));
        add_matrix_form(new MatrixFormHelmholtzEquation_imag_imag(1, 1, eps, mu, omega));

        add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_imag(0, 1, BDY_IMPEDANCE, beta));
        add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_real(1, 0, BDY_IMPEDANCE, beta));
    };

private:
    class MatrixFormHelmholtzEquation_real_real : public WeakForm::MatrixFormVol
    {
    public:
        MatrixFormHelmholtzEquation_real_real(int i, int j, double eps, double mu) :
            WeakForm::MatrixFormVol(i, j), eps(eps), mu(mu) {
        }

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }

    private:
        // Memebers.
        double eps;
        double mu;
    };

    class MatrixFormHelmholtzEquation_real_imag : public WeakForm::MatrixFormVol
    {

    private:
        // Memebers.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_real_imag(int i, int j, double mu, double omega, double sigma) :
            WeakForm::MatrixFormVol(i, j), mu(mu), omega(omega), sigma(sigma) {
        }

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return -omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
    };

    class MatrixFormHelmholtzEquation_imag_real : public WeakForm::MatrixFormVol
    {
    private:
        // Memebers.
        double mu;
        double omega;
        double sigma;

    public:
        MatrixFormHelmholtzEquation_imag_real(int i, int j, double mu, double omega, double sigma) :
            WeakForm::MatrixFormVol(i, j), mu(mu), omega(omega), sigma(sigma) {
        }

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return  omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }

    };

    class MatrixFormHelmholtzEquation_imag_imag : public WeakForm::MatrixFormVol
    {
    public:
        MatrixFormHelmholtzEquation_imag_imag(int i, int j, double eps, double mu, double omega) :
            WeakForm::MatrixFormVol(i, j), eps(eps), mu(mu), omega(omega) {
        }

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
        }
        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }

        // Memebers.
        double eps;
        double mu;
        double omega;
    };


    class MatrixFormSurfHelmholtz_real_imag : public WeakForm::MatrixFormSurf
    {

    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_real_imag(int i, int j, std::string area, double beta)
            : WeakForm::MatrixFormSurf(i, j, area), beta(beta){ }

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return beta * int_u_v<Real, Scalar>(n, wt, u, v);
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
    };



    class MatrixFormSurfHelmholtz_imag_real : public WeakForm::MatrixFormSurf
    {

    private:
        double beta;
    public:
        MatrixFormSurfHelmholtz_imag_real(int i, int j, std::string area, double beta)
            : WeakForm::MatrixFormSurf(i, j, area), beta(beta){ }

        template<typename Real, typename Scalar>
        Scalar matrix_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
        {
            return - beta*int_u_v<Real, Scalar>(n, wt, u, v);
        }
        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
        }
    };
};


