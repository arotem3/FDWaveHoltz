#include "WaveHoltz.hpp"
#include "linalg.hpp"

#include <armadillo>

using namespace wh;

static void initial_condition(double& u, double& v, double x, double w)
{
    const double sp = std::sin(M_PI * x);
    const double s2p = std::sin(2 * M_PI * x);
    const double sw = std::sin(w * x);
    const double cw = std::cos(w * x);
    u =  2.0 * sp * sp * sw;
    v = -2.0 * (w * cw * sp * sp + M_PI * s2p * sw);
}

static dvec linspace(double a, double b, int n)
{
    const double h = (b - a) / (n - 1);
    
    dvec x(n);
    for (int k=0; k < n; ++k)
    {
        x(k) = a + h * k;
    }

    return x;
}

struct factorization
{
    arma::cx_mat q;
    arma::cx_mat r;

    factorization(const arma::cx_mat& a)
    {
        arma::qr(q, r, a);
    }

    arma::cx_vec inv(const arma::vec& x)
    {
        return arma::solve(arma::trimatu(r), q.t() * x);
    }
};

int main()
{
    constexpr double alpha = (2.0 * M_PI * M_PI - 3.0) / (12.0 * M_PI);

    const double omega_start = 10, omega_end = 30, omega_delta = 0.5;
    const double a = 0.0, b = 2.0;

    std::cout << std::fixed << std::setprecision(4);

    std::ofstream out("solution/convergence_rate1d.txt");
    out << std::setprecision(10);
    out << "omega,e,mu,min parabolic distance,kappa\n";

    for (double w = omega_start; w <= omega_end; w += omega_delta)
    {
        const double omega = w * M_PI;
        const int n = num_points_per_unit_length(omega + 2*M_PI);

        const dvec x = linspace(a, b, n);
        const double h = x(1) - x(0);

        const char bc[] = "no";
        waveholtz1d WH(omega, &n, h, bc);

        arma::mat A(2*n, 2*n);
        WH.get_wave_op().as_matrix( A.memptr() );
        arma::cx_vec eigvals;
        arma::cx_mat R;
        arma::eig_gen(eigvals, R, A, "balance");
        
        const double kappa = arma::cond(R);

        eigvals /= omega;
        arma::vec parabolic_distance = -arma::real(eigvals) + alpha * arma::square(arma::imag(eigvals) - 1.0);
        const double min_distance = parabolic_distance.min();

        // pre-factor R
        factorization qr(R);

        dmat U(n, 2);
        for (int i=0; i < n; ++i)
        {
            initial_condition(U(i, 0), U(i, 1), x(i), omega);
        }

        arma::cx_vec mu = qr.inv(arma::vec(U, 2*n, false, false));

        const double e0 = norm(U);
        const double mu0 = arma::norm(mu);

        WH.S(U);
        mu = qr.inv(arma::vec(U, 2*n, false, false));

        const double e1 = norm(U);
        const double mu1 = arma::norm(mu);

        std::cout << "omega = " << w << " pi | n = " << std::setw(5) << n << " | min distance = " << min_distance << " | e1/e0 = " << std::setw(8) << e1/e0 << " | mu1/mu0 = " << std::setw(8) << mu1/mu0 << "\n";
        out << w << ", " << e1 / e0 << ", " << mu1 / mu0 << ", " << min_distance << ", " << kappa << "\n";
    }

    return 0;
}