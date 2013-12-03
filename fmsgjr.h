// fmsvalue.h - Valuation of options
// F = f exp(-kappa(t,s) + s Z_t) E exp(s Z_t) = exp(kappa(t,s))
// Model class must implement cdf(k) = P(F <= k) and cdf_(k) = E F/f 1(F <= k).
// E max{k - F, 0} = k P(F <= k) - f P(F/f 1( F <= k).
#pragma once
#pragma warning(disable: 4244)
#include "../xll8/xll/ensure.h"
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

	namespace black {
		// df .. dt get incremented and should be set to 0 on first use
		// return option value and update delta, gamma, vega, and theta
		template<class X>
		inline X greeks(X f, X s, X k, X t, X* df = 0, X* ddf = 0, X* ds = 0, X* dt = 0)
		{
			ensure (f >= 0);
			ensure (s >= 0);
			ensure (t >= 0);

			X c = 1;
			// negative strike means put
			if (k < 0) {
				c = -1;
				k = -k;
			}

			// boundary cases
			if (f == 0 || s == 0 || t == 0) {
				if (df) 
					*df += 1.*(c*f > c*k);
				if (ddf && f == k) 
					*ddf += std::numeric_limits<X>::infinity(); // really delta function at k

				return c*f > c*k ? c*(f - k) : 0;
			}

			if (k == 0) {
				if (df) 
					*df += 1;

				return f;
			}

			X srt = s*sqrt(t);
			X d2 = log(f/k)/srt - srt/2;
			X d1 = d2 + srt;
			X Nd1 = distribution::standard_normal<X>::cdf(c*d1);
			X Nd2 = distribution::standard_normal<X>::cdf(c*d2);
			X nd1(0);

			if (ddf || ds || dt)
				nd1 = distribution::standard_normal<X>::pdf(d1);

			if (df)
				*df += c*Nd1;

			if (ddf)
				*ddf += nd1/(f*srt);

			if (ds)
				*ds += f*srt*nd1/s;

			if (dt)
				*dt += -f*srt*nd1/(2*t); // negative of dv/dt

			return c*(f*Nd1 - k*Nd2);
		}

		// Black implied volatility from call value v (v > 0) or put value (v < 0).
		template<class X = double>
		inline X implied_volatility(const X& f, const X& v, const X& k, const X& t, X s0 = 0, const X& eps = 1e-8, size_t max_iter = 100)
		{
			// ensure price in 0 - infty vol range
			ensure (f > 0);
			ensure (v > (std::max)(f - k, 0.));
			ensure (v < f);
			ensure (k > 0);
			ensure (t > 0);

			if (s0 == 0) {
				// Brenner Subrahmanyam formula
				s0 = M_SQRT_2PI*v/(f*sqrt(t));
			}

			X ds(0);
			X v0 = greeks<X>(f, s0, k, t);
			if (fabs(v0 - v) < eps)
				return s0; //  got lucky

			// bracket the root
			X alpha = 2;
			if (v0 > v) {
				while (v0 > v) {
					s0 /= alpha;
					v0 = greeks<X>(f, s0, k, t);
				}
			}
			else {
				while (v0 < v) {
					s0 *= alpha;
					v0 = greeks<X>(f, s0, k, t);
				}
			}

			// Newton Raphson
			X dv = v0 - v;
			while (fabs(dv) > eps && max_iter--) {
				ds = 0;
				v0 = greeks<X>(f, s0, k, t, 0, 0, &ds);
				dv = (v0 - v);
				s0 = s0 - dv/ds;
			}

			return s0;
		}
	} // black

} // fms