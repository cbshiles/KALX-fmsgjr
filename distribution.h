// distribution.h
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "polynomial.h"

// pi
#ifndef M_PI
#define M_PI       3.141592653589793238462643383279502884197169399375105820974944L				   
#endif
// 2/pi
#ifndef M_2_PI
#define M_2_PI     0.636619772367581343075535053490057448137838582961825794990669L
#endif
// sqrt(2*pi)
#ifndef M_SQRT_2PI
#define M_SQRT_2PI 2.506628274631000502415765284811045253006986740609938316629923L
#endif
// sqrt(2)
#ifndef M_SQRT_2
#define M_SQRT_2    1.414213562373095048801688724209698078569671875376948073176679L
#endif

namespace fms {

	namespace distribution {

		// cumulants of the Esscher transform
		template<class X = double>
		void esscher(size_t nk, const X* kappa, const X& gamma, X* kappa_)
		{
			ensure (kappa + nk <= kappa_ || kappa_ + nk < kappa);

			if (nk > 0)
				std::copy(kappa, kappa + nk, kappa_);

			for (size_t k = 0; k < nk; ++k) {
				X gn = gamma; // gamma^n/n!
				for (size_t n = 0; n + k < nk; ++n, gn *= gamma/n)
					kappa_[k] += kappa[n + k]*gn;
			}
		}

		template<class X = double>
		struct standard_normal {
			static X pdf(const X& x)
			{
				return exp(-x*x/2)/M_SQRT_2PI;
			}
			// n-th derivative
			static X cdf(const X& x, size_t n = 0)
			{
				if (n == 0) 
					return (1 + erf(x/M_SQRT_2))/2;

				// n&~1 = n%2 ???
				return (n%2 == 1 ? 1 : -1)*polynomial::Hermite(n - 1, x)*exp(-x*x/2)/M_SQRT_2PI;
			}
			template<class P = double>
			static P inv(P p, P x = std::numeric_limits<P>::quiet_NaN(), P eps = sqrt(std::numeric_limits<P>::epsilon()))
			{
				if (_isnan(x))
					x = tan(M_PI*(p - .5))*sqrt(M_2_PI)*pow(4*p*(1-p), M_2_PI);

				//!!! if (p < 0.05 || p > 0.95)

				for (P x_ = x + 1; fabs(x_ - x) > eps; std::swap(x_,x)) {
					x_ = x - (cdf(x) - p)/pdf(x);
				}

				return x;
			}

			// alias
			static X F(const X& x, size_t n = 0)
			{
				return cdf(x, n);
			}
			// perturb first k cumulants by kappa
			template<class X = double, class K>
			static X F_(const X& z, size_t n = 0, size_t k = 0, const K* kappa = 0, size_t N = 0, X* dG = 0)
			{
				X g = F(z, n);

				if (k == 0)
					return g;

				X _1(-1); // (-1)^n
				// put in policy class
				X tol(0);
				X eps(1e-8);
				for (size_t i = 1; i < 30 && (tol == 0 || fabs(tol) > eps); ++i, _1 *= -1) {
					X Bi = Bell(i, k, kappa, true);

					if (Bi) {
						X dg = _1*Bi*F(z, i + n);
						if (i < N)
							dG[i] = dg;
						g += dg;

						if (tol && dg)
							tol = (std::min)(tol, fabs(dg));
						else if (!tol && dg)
							tol = fabs(dg);
					}
				}

				return g;
			}
		};
		
	} // distribution

}