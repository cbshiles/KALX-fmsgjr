// fms_distribution.h - probability distributions
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

// sqrt(2*pi)
#define M_SQRT_2PI 2.506628274631000502415765284811045253006986740609938316629923
// log_2(e)
#define M_LOG2_E   1.442695040888963407359924681001892137426645954152985934135449

namespace fms {
	namespace distribution {

		// probability distribution function, P(X = x) dx
		template<class D, class T = double>
		inline T pdf(D d, T x);
		template<class I, class D, class T = double>
		inline T cdf(D d, T x);
		template<class I, class D, class T = double>
		inline T inv(D d, T x);

		template<class M = double, class S = double>
		struct normal {
			typedef M mean_type;
			typedef S standard_deviation_type;
			M mean;
			S standard_deviation;
			normal(M m = 0, S s = 1)
				: mean(m), standard_deviation(s)
			{ }
		};

		template<class T = double, class M = double, class S = double>
		inline T pdf(normal<M,S> d, T t)
		{
			decltype(t - d.mean) m = t - d.mean;
			S s = d.standard_deviation;

			return exp(-m*m/(2*s*s))/(s*M_SQRT_2PI);
		}

	}

}