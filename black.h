// black.h - Black forward option value.
#pragma once
#include "distribution.h"
#include "../fmsroot/root1d.h"

namespace fms {
	
	namespace black {

		// P(cF > ck)
		template<class X>
		inline X binary(X f, X s, X k, X t)
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
				return c*f > c*k ? 1 : 0;
			}

			if (k == 0) {
				return c == 1 ? 1 : 0;
			}

			X srt = s*sqrt(t);
			X d2 = log(f/k)/srt - srt/2;
			X Nd2 = distribution::standard_normal<X>::cdf(c*d2);

			return Nd2;
		}

		// E max{c(F - k), 0}
		template<class X>
		inline X value(X f, X s, X k, X t)
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
				return c*f > c*k ? c*(f - k) : 0;
			}

			if (k == 0) {
				return c == 1 ? f : 0;
			}

			X srt = s*sqrt(t);
			X d2 = log(f/k)/srt - srt/2;
			X d1 = d2 + srt;
			X Nd1 = distribution::standard_normal<X>::cdf(c*d1);
			X Nd2 = distribution::standard_normal<X>::cdf(c*d2);

			return c*(f*Nd1 - k*Nd2);
		}

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

		template<class X>
		inline X corrado_miller_implied_volatility(const X& f, const X& v, const X& k, const X& t)
		{
			X s = v - (f - k)/2;
			s *= s;
			s -= (f - k)*(f - k)/M_PI;
			s = sqrt(s);
			s += v - (f - k)/2;
			s *= M_SQRT_2PI/((f + k)*sqrt(t));

			return s;
		}

		// Black implied volatility from call value v.
		template<class X = double>
		inline X implied_volatility(const X& f, X v, X k, const X& t, X s0 = 0)
		{
			if (k < 0) {
				// put to call
				k = -k;
				//c - p = f - k;
				//c = p + f - k;
				v = v + f - k;
			}
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
			
			// vega must be at least 1
			X ds(0);
			black::greeks<X>(f, s0, k, t, 0, 0, &ds);
			while (ds < 1) {
				s0 *= 1.6;
				ds = 0;
				black::greeks<X>(f, s0, k, t, 0, 0, &ds);
			}

			fms::root1d<X> r([f, v, k, t](X s) { return black::value<X>(f, s, k, t) - v; },
							 [f, v, k, t](X s) { X ds(0); black::greeks<X>(f, s, k, t, 0, 0, &ds); return ds; }
			);
			r.push(s0);
			while (!r.residual(1e-8) || !r.delta(1e-8)) {
				r.newton();
			}			

			return r.x[0];
		}
	} // black
} // fms