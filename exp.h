// exp.h - exponential function
#pragma once

namespace fms {
	namespace math {
		template<class T = double>
		inline T exp(T x)
		{
			double n = 1;
			T e = 0, e_ = 1, xn = 1;

			while(e_ != e)
				e_ = (e = e_) + (xn*=(x/n++));

			return e_;
		}
	}
}