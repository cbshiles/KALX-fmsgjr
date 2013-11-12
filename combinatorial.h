// combinatorial.h
#pragma once
#ifndef ensure
#include <cassert>
#define ensure assert
#endif
#include <map>
#include <vector>

namespace fms {

	namespace combinatorial {
		// memoized n!
		inline long long factorial(size_t n)
		{
			static std::vector<long long> n_ = {
				1, 1, 2, 6, 24, 120, 
				720, 5040, 40320, 362880, 3628800, 
				39916800, 479001600, 6227020800, 87178291200, 1307674368000, 
				20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000
			};
			
			if (n >= n_.size()) {
				size_t j = n_.size();
				n_.resize(n + 1, 0);
				for (size_t i = j; i <= n; ++i)
					n_[i] = i*n_[i-1];
			}

			return n_[n];
		}

		// memoized C(n,k) = n!/k!(n - k)!
		inline size_t choose(size_t n, size_t k) {
			ensure (k <= n);
			static std::map<std::pair<size_t,size_t>,size_t> cnk_; // = {{{0,0},0},...

			if (k == 0 || k == n) {
				return 1;
			}
			else {
				auto nk = std::make_pair(n,k);
				auto i = cnk_.find(nk);

				if (i != cnk_.end())
					return i->second;
				
				return cnk_[nk] = choose(n - 1, k - 1) + choose(n - 1, k);
			}
		};

	} // combinatorial
} // fms