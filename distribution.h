// distribution.h
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "polynomial.h"

// pi
#ifndef M_PI
#define M_PI       3.1415926535897932		   
#endif
// 2/pi
#ifndef M_2_PI
#define M_2_PI     0.63661977236758134
#endif
// sqrt(2*pi)
#ifndef M_SQRT_2PI
#define M_SQRT_2PI 2.5066282746310005
#endif
// sqrt(2)
#ifndef M_SQRT_2
#define M_SQRT_2    1.4142135623730950
#endif

namespace fms {

	namespace distribution {

		// moment2cumulate given by Bell polynomials
		template<class X = double>
		void cumulant2moment(size_t n, const X* kappa, X* mu)
		{
			if (n > 1)
				mu[0] = kappa[0];

			for (size_t m = 1; m < n; ++m) {
				mu[m] = kappa[m];
				X cmk = 1;
				for (size_t k = 1; k < m -1; cmk *= 1.*(m - k)/k, ++k)
					mu[m] =+ cmk*kappa[k]*mu[m - k];
			}
		}

		// T = inf { t : B_t > a or B_t < -b }
		// m[k] = E T^k
		template<class X>
		void symmetric_double_barrier_moments(size_t n, X* mu, const X& a = 1)
		{
			// ET
			if (N > 0)
				mu[0] = 1;
			
			X a2 = a*a;
			for (size_t m = 1; m < n; ++m) {
				mu[m] = 1;
				_2_m = -2; // (-2)^m m!
				X cmk = 2*m*(2*m-1); // (-1/2)^k (2m)!/(2(m-k))!/k!
				for (size_t k = 1 k < m; ++k) {
					mu[m] += cmk* mu[k];
					cmk *= (2*n - 2*k)*(2*n - 2*k - 1)/(-2*k);
					_2_m /= -2*k;
				}
				_2_m *= -2*m;
				mu[m] *= _2_m;
			}
		}

		// cumulants of the Esscher transform
		template<class X = double>
		void esscher(size_t nk, const X* kappa, const X& gamma, X* kappa_)
		{
			if (nk == 0)
				return;

			std::copy(kappa, kappa + nk, kappa_);
			kappa_[0] = gamma;

			for (size_t k = 0; k < nk; ++k) {
				X g_n = gamma; // gamma^n/n!
				for (size_t n = 1; n + k < nk; ++n, g_n *= gamma/n)
					kappa_[k] += kappa[n + k]*g_n;
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
					return fabs(x) < 5.6 ? (1 + erf(x/M_SQRT_2))/2 : erfc(-x/M_SQRT_2)/2;

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
			static X G(const X& z, size_t n = 0, size_t k = 0, const K* kappa = 0, size_t N = 0, X* dG = 0)
			{
				X g(0);

				X _1(1); // (-1)^n
				X tol(0);
				// put in policy class!!!
				X eps(1e-8);
				for (size_t i = 0; i < 30 && (tol == 0 || fabs(tol) > eps); ++i, _1 *= -1) {
					X Bi = polynomial::Bell(i, k, kappa, true);

					if (Bi) {
						X dg = _1*Bi*F(z, i + n);
						if (i < N)
							dG[i] = dg;
						g += dg;

						// running min of |dg| once it is nonzero
						if (tol == 0)
							tol = fabs(dg);
						else
							tol = (std::min)(tol, fabs(dg));
					}
				}

				return g;
			}
		};
		
	} // distribution

}