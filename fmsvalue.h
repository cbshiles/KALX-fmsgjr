// fmsvalue.h - Valuation of options
// F = f exp(-kappa(t,s) + s Z_t) E exp(s Z_t) = exp(kappa(t,s))
// Model class must implement cdf(k) = P(F <= k) and cdf_(k) = E F/f 1(F <= k).
// E max{k - F, 0} = k P(F <= k) - f P(F/f 1( F <= k).
#pragma once
#pragma warning(disable: 4244)
#ifndef ensure
#include <cassert>
#define ensure assert
#endif
#include <algorithm>
#include <functional>
#include <map>
#include <utility>
#include <valarray>
#include <vector>
//#include "combinatorial.h"
#include "distribution.h"
#include "instrument.h"
#include "model.h"
#include "polynomial.h"

namespace fms {

	template<class F, class S, class K, class T, template<class,class> class M, template<class,class> class I>
	inline auto value(const M<F,S>& m, const I<K,T>& o)
		-> decltype(m.forward + m.volatility + o.strike + o.expiration)
	{
		F f = m.forward;
		S s = m.volatility;
		K k = o.strike;
		T t = o.expiration;

		// boundary cases
		if (1 + f == 1) // zero underlying
			return k > 0 ? 0 : -k;
		if (1 + s == 1 || 1 + t == 1) // intrinsic
			return k > 0 ? (std::max)(f - k, 0.) : (std::max)(-k - f, 0.);
		if (1 + k == 1)
			return k > 0 ? f : 0;

		return k > 0 ? f*m.cdf(k,t) - k*m.cdf_(k,t) : -k*m.cdf(-k,t) - f*m.cdf_(-k,t);
	}

} // namespace fms