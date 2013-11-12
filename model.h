// model.h
#pragma once
#include <valarray>
#include "distribution.h"

namespace fms {
	namespace model {

 		// f exp(-s^2t/2 + sB_t)
 		template<class F = double, class S = double>
		struct black {
			F forward;
			S volatility;
			black(const F& f, const S& s)
				: forward(f), volatility(s)
			{ }
			black(const black& b) = default;
			black& operator=(const black& b) = default;
			~black()
			{ }

			// P(F <= k) = E[1(F <= k)]
			template<class K = double, class T = double>
			auto cdf(K k, T t) const -> decltype(forward + volatility + k + t)
			{
				F f = forward;
				S s = volatility;
				auto gamma = s*sqrt(t);
				auto z = (log(k/f) + gamma*gamma/2)/gamma;

				return distribution::standard_normal<decltype(z)>::cdf(z);
			}
			// E[F 1(F <= k)]/E[F]
			template<class K = double, class T = double>
			auto cdf_(K k, T t) const -> decltype(forward + volatility + k + t)
			{
				F f = forward;
				S s = volatility;
				auto gamma = s*sqrt(t);
				auto z = (log(k/f) - gamma*gamma/2)/gamma;

				return distribution::standard_normal<decltype(z)>::cdf(z);
			}
		};

		// Generalized Jarrow-Rudd model
		// !!!use iterator for cumulants instead of array!!!
		// f exp(-s^2t/2 + sX_t) kappa(X_i) are cumulant
 		template<class F = double, class S = double>
		class gjr {
			bool black_; // all cumulants are 0
		public:
			F forward;
			S volatility;
			std::valarray<F> kappa;
			mutable std::valarray<F> kappa_;

			gjr(const F& f = 0, const S& s = 0, size_t k = 0, const double* c = 0)
				: forward(f), volatility(s), kappa(c, k), kappa_(k)
			{
				black_ = k == 0 || ((kappa.max)() == 0 && (kappa.min)() == 0);

				at(1);
			}
			gjr(const gjr&) = default;
			gjr& operator=(const gjr&) = default;
			~gjr()
			{ }

			// set X_t
			template<class T>
			void at(const T& t) const
			{
				auto gamma = volatility*sqrt(t);

				distribution::esscher(kappa.size(), &kappa[0], gamma, &kappa_[0]);
			}

			// P(F <= k) = E[1(F <= k)]
			template<class K = double, class T = double>
			auto cdf(K k, T t) const -> decltype(forward + volatility + k + t)
			{
				auto f = forward;
				auto s = volatility;
				auto gamma = s*sqrt(t);
				auto z = (log(k/f) + gamma*gamma/2)/gamma;

				return distribution::standard_normal<decltype(z)>::F_(z, 0, kappa.size(), &kappa[0]);
			}

			// P(F <= k) = E[1(F <= k)]
			template<class K = double, class T = double>
			auto cdf_(K k, T t) const -> decltype(forward + volatility + k + t)
			{
				auto f = forward;
				auto s = volatility;
				auto gamma = s*sqrt(t);
				auto z = (log(k/f) - gamma*gamma/2)/gamma;

				at(t);

				return distribution::standard_normal<decltype(z)>::F_(z, 0, kappa_.size(), &kappa_[0]);
			}
			
		};

	} // model
} // fms