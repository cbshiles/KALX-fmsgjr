// fmsvalue.h - Valuation of options
// F = f exp(-kappa(t,s) + s Z_t) E exp(s Z_t) = exp(kappa(t,s))
// Model class must implement cdf(k) = P(F <= k) and cdf_(k) = E F/f 1(F <= k).
// E max{k - F, 0} = k P(F <= k) - f P(F/f 1( F <= k).
#pragma once
#pragma warning(disable: 4244)
#include <cassert>
#define ensure(x) assert(x)
//#include "../xll8/xll/ensure.h"
#include <algorithm>
#include <functional>
#include <map>
#include <utility>
#include <valarray>
#include <vector>
#include "combinatorial.h"
#include "cumulant.h"
#include "distribution.h"
#include "instrument.h"
#include "model.h"
#include "polynomial.h"

namespace fms {

	template<class F, class S, class K, class T, template<class,class> class M>
	inline auto put(const F& f, const S& s, const K& k, const T& t, const M<F,S>& m = model::black<F,S>(f,s))
		-> decltype(f + s + k + t)
	{
		// boundary cases
		if (1 + f == 1) // zero underlying
			return k > 0 ? 0 : -k;
		if (1 + s == 1 || 1 + t == 1) // intrinsic
			return (std::max)(-k - f, 0.);
		if (1 + k == 1)
			return 0;

		return k < 0 ? -k*m.cdf(-k,t) - f*m.cdf_(-k,t) : k*m.cdf(k,t) - f*m.cdf_(k,t);
	}
	template<class F, class S, class K, class T, template<class,class> class M>
	inline auto call(const F& f, const S& s, const K& k, const T& t, const M<F,S>& m = model::black<F,S>(f,s))
		-> decltype(f + s + k + t)
	{
		return f - k + put(f, s, k, t, m);
	}
	template<class F, class S, class K, class T, template<class,class> class M>
	inline auto value(const F& f, const S& s, const K& k, const T& t, const M<F,S>& m = model::black<F,S>(f,s))
		-> decltype(f + s + k + t)
	{
		return k < 0 ? put(f, s, k, t, m) : call(f, s, k, t, m);
	}

	template<class F, class S, class K, class T, template<class,class> class M, template<class,class> class I>
	inline auto put(const M<F,S>& m, const I<K,T>& o)
		-> decltype(m.forward + m.volatility + o.strike + o.expiration)
	{
		return put(m.forward, m.volatility, o.strike, o.expiration, m);
	}
	template<class F, class S, class K, class T, template<class,class> class M, template<class,class> class I>
	inline auto call(const M<F,S>& m, const I<K,T>& o)
		-> decltype(m.forward + m.volatility + o.strike + o.expiration)
	{
		return call(m.forward, m.volatility, o.strike, o.expiration, m);
	}
	template<class F, class S, class K, class T, template<class,class> class M, template<class,class> class I>
	inline auto value(const M<F,S>& m, const I<K,T>& o)
		-> decltype(m.forward + m.volatility + o.strike + o.expiration)
	{
		return o.strike < 0 ? put(m, o) : call(m, o);
	}

} // fms