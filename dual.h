// dual.h - dual numbers
#pragma once
#include <initializer_list>
#include <valarray>

namespace dual {

	// dual number with maximum size N
	template<class X = double>
	class number {
		size_t N_;
		std::valarray<X> x_; // x[0] + x[1] J + x[2] J^2 + ...
	public:
		number(size_t N = 0)
			: N_(N), x_{0}
		{ }
		number(const X& x, size_t N = 0)
			: N_(N), x_(x, 1)
		{ }
		number(std::initializer_list<X> x, size_t N = 0)
			: N_(N), x_(x)
		{ }
		number(const number&) = default;
		number& operator=(const number&) = default;
		number(number&& x)
			: N_(x.N_), x_(std::move(x))
		{ }
		number& operator=(number&& x)
		{
			if (this != & x) {
				N_ = x.N_;
				x_ = std::move(x);
			}

			return *this;
		}
		~number()
		{ }

		operator std::valarray<X>&()
		{
			return x_;
		}
		operator const std::valarray<X>&() const
		{
			return x_;
		}

		size_t size(void) const
		{
			return x_.size();
		}
		size_t order(void) const
		{
			return N_;
		}

		const X& operator[](size_t i) const
		{
			return x_[i];
		}
		X& operator[](size_t i)
		{
			return x_[i];
		}
		// f(x) = f(0)I + f'(0)J + f''(0)J^2/2! + ...
		X operator()(size_t i) const
		{
			return x_[i]*factorial(i);
		}

		number& operator+=(const X& x)
		{
			x_[0] += x;

			return *this;
		}
		number& operator+=(const number& y)
		{
			if (size() < y.size()) {
				resize(y.size());
			}

			x_[std::slice(0, y.size(), 1)] += y;

			return *this;
		}
		// x y = x[0] y[0] I + (x[0] y[1] + x[1] y[0]) J + ...
		number& operator*=(const number& y)
		{
			size_t n = size() + y.size() - 1; // natural size

			if (order() || y.order())
				n = std::max(order(), y.order());

			std::valarray<X> z(X(0), n);

			for (size_t i = 0; i < z.size(); ++i) {
				for (size_t j = 0; j < x_.size(); ++j) {
					if (i - j < y.size())
						z[i] += x_[j]*y[i - j];
				}
			}

			z.swap(x_);

			return *this;
		}
		number& operator/=(const X& y)
		{
			x_ /= y;

			return *this;
		}
		// (a + bJ)^-1 = 1/a - b/a^2 J + b^2/a^3 J^2 + ...
		number inv(void) const
		{
			X a(x_[0]);
			number y(1/a, N_);
			number b(*this);

			b.shift(1);
			number c(1/a, N_);
			for (size_t i = 1; i < N_; ++i) {
				c.shift(-1);
				c *= b;
				c /= -a;
				y += c;
			}

			return y;
		}
	protected:
		number& resize(size_t n)
		{
			if (N_)
				n = (std::min)(n, N_);

			if (n != size()) {
				std::valarray<X> x(x_);
				x_.resize(n);
				if (n > x.size())
					x_[std::slice(0, x.size(), 1)] = x;
				else
					x_ = x[std::slice(0, n, 1)];
			}

			return *this;
		}
		number& shift(int n)
		{
			if (n > 0) {
				x_ = x_.shift(n);
				// ??? resize ???
			}
			else if (n < 0) {
				std::valarray<X> x(x_.size() - n);
				x[std::slice(-n, x_.size() - n, 1)] = x_;
				x_.swap(x);
			}

			return *this;
		}
		private:
			static double factorial(size_t n)
			{
				static double n_[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 393628800,
					39916800, 479001600, 6227020800, 87178291200, 1307674368000,
					20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 
					2432902008176640000};

				ensure (n < sizeof(n_)/sizeof(*n_));

				return n_[n];
			}
	};

	template<class X>
	inline number<X>op_plus() {}

} // dual

template<class X>
inline dual::number<X> operator+(const dual::number<X>& x, const X& y)
{
	dual::number<X> z(x);

	return z += y;
}
template<class X>
inline dual::number<X> operator+(const X& x, const dual::number<X>& y)
{
	dual::number<X> z(y);

	return z += x;
}
template<class X>
inline dual::number<X> operator+(const dual::number<X>& x, const dual::number<X>& y)
{
	dual::number<X> z(x);

	return z += y;
}

template<class X>
inline dual::number<X> operator*(const dual::number<X>& x, const dual::number<X>& y)
{
	dual::number<X> z(x);

	return z *= y;
}