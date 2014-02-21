// combinatorial.h
#pragma once
#ifndef ensure
#include <cassert>
#define ensure assert
#endif
#include <iterator>
#include <limits>
#include <map>
#include <vector>

namespace fms {

	namespace combinatorial {

		class factorial_iterator : public std::iterator<std::forward_iterator_tag, long long> {
			long long i_, n_;
		public:
			factorial_iterator()
				: i_(0), n_(1)
			{ }
			~factorial_iterator()
			{ }

			long long operator*() const
			{
				return n_;
			}
			factorial_iterator& operator++()
			{
				++i_;
				n_ *= i_;

				return *this;
			}
			factorial_iterator operator++(int)
			{
				factorial_iterator fi(*this);

				++i_;
				n_ *= i_;

				return fi;
			}


		};

		// memoized n!
		inline long long factorial(size_t n)
		{
			static std::vector<unsigned long long> n_ = {
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
			static std::vector<std::vector<unsigned long long>> cnk_ = {
				{1},
				{1, 1},
				{1, 2, 1}
			};

			if (n >= cnk_.size()) {
				size_t n_ = cnk_.size();
				cnk_.resize(n + 1);
				for (size_t m = n_; m <= n; ++m) {
					cnk_[m].resize(m + 1);
					cnk_[m][0] = 1;
					cnk_[m][m] = 1;
					for (size_t j = 1; j < m; ++j) {
						unsigned long long cmj_ = cnk_[m-1][j-1];
						unsigned long long cmj = cnk_[m-1][j];
						ensure (cmj_ < (std::numeric_limits<unsigned long long>::max)() - cmj);
						cnk_[m][j] = cmj_ + cmj;
					}
				}
			}
				
			return cnk_[n][k];
		};

	} // combinatorial
} // fms