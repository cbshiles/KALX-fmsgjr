// daly.h - Daly's improvements to Marsaglia's routine
// From: bill.daly@tradition-ny.com (Bill Daly)
// Message-ID: <4fa0d5e.0402241117.53c8da2a@posting.google.com>
// Message-ID: <4fa0d5e.0402181449.6dd203db@posting.google.com>
// standard normal cumulative distribution function
// cdf(normal<daly>, x)
#pragma once
#include <cmath>

namespace fms {
	namespace distribution {

		class daly {};

		template<class T = double, class M = double, class S = double>
		inline T cdf(normal<M,S>, T x)
		{
			I i = 1;
			T s=x,t=0,b=x,x2=x*x;

			while(s!=t)
				s=(t=s)+(b*=x2/(i+=2));

			return  static_cast<T>(.5+s*exp(-.5*x2-.91893853320467274178L));
		}

	} // distribution
} // fms