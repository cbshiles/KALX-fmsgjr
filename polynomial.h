// polynomial.h
#pragma once

namespace fms{
	namespace polynomial {



		template<class X>
		inline X Bell(size_t n, size_t nx, const X* x, bool reduced = false)
		{
			std::vector<X> B(n + 1, 0);

			B[0] = 1;

			for (size_t k = 1; k <= n; ++k) {
				X ckj = reduced ? 1./k : 1; // (k - 1 choose j)
				for (size_t j = 0;j < k && j < nx; ++j, ckj /= j) {
					B[k] += ckj * B[k - 1 - j]*x[j];
					if (!reduced) {
						ckj *= (k - j - 1);
					}
				}
			}

			return B[n];
		}

		// H_0(x) = 1, H_1(x) = x, H_{n+1}(x) = xH_n(x) - n H_{n-1}(x)
		template<class X>
		inline X Hermite(size_t n, const X& x)
		{
			std::vector<X> H(n+1, 0);

			if (n <= 0)
				return 1;
			if (n <= 1)
				return x;

			H[0] = 1;
			H[1] = x;

			for (size_t i = 2; i <= n; ++i)
				H[i] = x*H[i-1] - (i - 1)*H[i-2];
			
			return H[n];
		}

	} // polynomial
} // fms