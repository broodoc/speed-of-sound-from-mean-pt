//Copyright 2026 Caleb Broodo
#include "Models.h"

double Models::operator()(double* x, double* par) {
    const double xx = x[0];

    switch (type_f) {

        case Function::LevyTsallis:
            // m - par[0], C - par[1], n - par[2], dN/dy - par[3]
            return par[3] * (par[2] - 1) * (par[2] - 2) /
                   (2 * 3.14159 * par[1] * par[2] *
                    (par[1] * par[2] + par[0] * (par[2] - 2))) *
                   std::pow(1 + (std::sqrt(xx * xx + par[0] * par[0]) - par[0]) /
                                   (par[1] * par[2]),
                            -par[2]) *
                   xx;

        case Function::Tsallis:
            // V - par[0], q - par[1], T - par[2]
            return (2.0 * std::pow(xx, 2.0) * par[0] / 39.478418) *
                   (std::pow((1.0 + (par[1] - 1) *
                                    (std::sqrt(std::pow(0.139570, 2.0) + std::pow(xx, 2.0)) / par[2])),
                             (-1.0 * par[1] / (par[1] - 1.0))) +
                    std::pow((1.0 + (par[1] - 1) *
                                    (std::sqrt(std::pow(0.493677, 2.0) + std::pow(xx, 2.0)) / par[2])),
                             (-1.0 * par[1] / (par[1] - 1.0))) +
                    2.0 * std::pow((1.0 + (par[1] - 1) *
                                          (std::sqrt(std::pow(0.938272, 2.0) + std::pow(xx, 2.0)) / par[2])),
                                   (-1.0 * par[1] / (par[1] - 1.0))));

        case Function::Hagedorn:
            // m - par[0], pT0 - par[1], b - par[2], constant - par[3]
            return (xx * xx / std::sqrt(par[0] * par[0] + xx * xx)) *
                   std::pow((1 + xx / par[1]), -par[2]) * par[3];

        case Function::TBW: {
            const double A      = par[0];
            const double m      = par[1];
            const double T      = par[2];
            const double q      = par[3];
            const double beta_s = par[4];
            const double pT     = x[0];
            if (T <= 0.0 || q <= 1.0 || beta_s < 0.0 || beta_s >= 1.0) return 0.0;

            const int n_r   = 4;
            const int n_phi = 5;
            const int n_y   = 5;
            const double R  = 1.0;
            const double dr   = R / n_r;
            const double dphi = 2.0 * TMath::Pi() / n_phi;
            const double Ywin = 0.10;
            const double dy   = (2.0 * Ywin) / n_y;

            const double mT = std::sqrt(pT * pT + m * m);
            double sum = 0.0;

            for (int ir = 0; ir < n_r; ++ir) {
                const double r  = (ir + 0.5) * dr;
                const double br = beta_s * (r / R);
                if (br >= 0.999999) continue;
                const double rho = std::atanh(br);

                for (int ip = 0; ip < n_phi; ++ip) {
                    const double phi = -TMath::Pi() + (ip + 0.5) * dphi;
                    const double cosphi = std::cos(phi);

                    for (int iy = 0; iy < n_y; ++iy) {
                        const double y  = -Ywin + (iy + 0.5) * dy;
                        const double cy = std::cosh(y);

                        const double expo = (mT * cy * std::cosh(rho)
                                           - pT * std::sinh(rho) * cosphi) / T;
                        const double base = 1.0 + (q - 1.0) * expo;
                        if (base <= 0.0) continue;

                        const double weight = std::pow(base, -1.0 / (q - 1.0));
                        sum += r * cy * weight;
                    }
                }
            }

            const double shape = mT * sum * dr * dphi * dy;
            return A * shape * pT;
        }

        case Function::FokkerPlanck: {
            // double m = par[0] (unused), A=par[1], b=par[2], c=par[3], d=par[4], T=par[5]
            const double A = par[1], b = par[2], c = par[3], d = par[4], T = par[5];
            const double arg1 = xx + b;
            const double arg2 = -d * xx / T;
            if (arg1 <= 0 || T <= 0 || !std::isfinite(arg1) || !std::isfinite(arg2))
                return 1e-12;
            const double result = A * std::pow(arg1, -c) * std::exp(arg2);
            return (std::isfinite(result) && result >= 0) ? result : 1e-12;
        }

        case Function::PowerLaw:
            // par[0]=cs2 exponent, par[1]=scale (your form)
            return par[1] * std::pow(xx / par[1], par[0]);

        case Function::ExpDecay:
            return xx * (par[0] + par[1] * xx * xx) / std::exp(par[2] * xx);

        case Function::RightKneeFunction: {
            const double n     = x[0];
            const double cs2   = par[0];
            const double knee  = par[1];
            const double sigma = par[2];

            const double norm = sigma * TMath::Sqrt(2.0 / TMath::Pi());
            const double exponent = TMath::Exp(-TMath::Power(n - knee, 2) /
                                               (2.0 * TMath::Power(sigma, 2)));
            const double erfc = TMath::Erfc((n - knee) / (TMath::Sqrt(2.0) * sigma));
            if (erfc < 1e-10) return 1e9;

            const double f = n - (norm * exponent / erfc);
            if (f <= 1e-5) return 1e9;

            return std::exp(cs2 * std::log(n / f));
        }

        case Function::RightKneeFunctionApproximation:
            return TMath::Power((xx / par[0]), par[1]);

        default:
            return 0.0;
    }
}

