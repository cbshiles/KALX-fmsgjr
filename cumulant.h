// cumulant.h - cumulants
#pragma once

namespace fms {
	namespace cumulant {

		template<class X = double>
		void scale(const X& c, size_t n, X* kappa)
		{
			X cn(c);

			for (size_t i = 0; i < n; ++i) {
				kappa[i] *= cn;
				cn *= c;
			}
		}

		// poisson with mean mu
		template<class X = double>
		void poisson(const X& mu, size_t n, X* kappa)
		{
			ensure (mu > 0);

			for (size_t i = 0; i < n; ++i) {
				kappa[i] = mu;
			}
		}

		// exponential with mean mu
		template<class X = double>
		void exponential(const X& mu, size_t n, X* kappa)
		{
			ensure (mu > 0);
			double mun(mu); // mu^n
			double n_(1); // (n - 1)!

			for (size_t i = 0; i < n; ++i, n_ *= i) {
				kappa[i] = n_*mun; // (n - 1)! mu^n
				mun *= mu;
			}
		}

		// gamma with mean mu = a/b variance s2 = a/b^2
		template<class X = double>
		void gamma(const X& mu, const X& s2, size_t n, X* kappa)
		{
			ensure (mu > 0);
			ensure (s2 > 0);

			X b(mu/s2);
			X a(mu*b);
			X bn(b);
			double n_(1); // (n - 1)!

			for (size_t i = 0; i < n; ++i, n *= i) {
				kappa[i] = n_*a/bn; // (n - 1)! alpha/beta^n
				bn *= b;
			}
		}

	} // cumulant
} // fms