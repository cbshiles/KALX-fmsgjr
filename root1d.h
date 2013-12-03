// root1d.h - One dimensional root finding
#pragma once
#include <functional>
#include <limits>
#include <deque>

namespace fms {
	namespace root1d {

		/*
		template<class X = double> struct state;

		template<class X = double>
		struct policy {
			std::function<void(void)> initialize;
			std::function<void(state&) update;
			std::function<bool(void)> solved;
		};
		*/
		
		template<class X = double>
		struct state {
			std::deque<X> x, y; // successive approximations
			X abs, rel; // absolute and relative error
			size_t iter; // maximum number of iterations
			state(X a = 0, X r = 0, size_t m = (std::numeric_limits<size_t>::max)())
				: abs(a), rel(r), iter(m)
			{}
			state(const state&) = default;
			state& operator=(const state&) = default;
			~state()
			{ }

			/*
			void initialize(...)
			void update()
			bool solved();
			*/

			//
			// convergence conditions
			//

			// use for roots close to 0
			bool test_inverval() const
			{
				X min_ = x[0]*x[1] < 0 ? 0 : (std::min)(fabs(x[0]), fabs(x[1]));

				return fabs(x[0] - x[1]) < abs + rel*min_;
			}
			bool test_delta() const
			{
				return fabs(x[0] - x[1]) < abs + rel*x[0];
			}
			bool test_residual() const
			{
				return fabs(y[0]) < abs;
			}
			bool bracketed() const
			{
				ensure (x.size() >= 2);

				return y[0]*y[1] < 0;
			}

			//
			// iteration methods
			//

			// assumes monotonic function
			X bracket(const std::function<X(X)>& f, X x0, X alpha = 2)
			{
				x.push_front(x0);
				y.push_front(f(x[0]));

				x.push_front(x[0]*alpha);
				y.push_front(f(x[0]));

				if (bracketed()) {					
					return x[0];
				}

				if (fabs(y[0]) > fabs(y[1])) {
					alpha = 1/alpha;

					x[0] = x[0]*alpha;
					y[0] = f(x[0]);

					ensure (fabs(y[0]) < fabs(y[1])); // monotonic
				}

				while (!bracketed()) {
					x.push_front(x[0]*alpha);
					y.push_front(f(x[0]));
				}

				return x[0];
			}
			X bisect(const std::function<X(X)>& f)
			{
				ensure (x.size() >= 2);
				ensure (x[0]*x[1] < 0);

				x.push_front((x[0] + x[1])/2);
				y.push_front(f([0]));

				if (y[0]*y[1] > 0) {
					std::swap(x[1], x[2]);
					std::swap(y[1], y[2]);
				}

				return x[0];
			}
			X secant(const std::function<X(X)>& f)
			{
				ensure (x.size() >= 2);
				ensure (y[0] != y[1]);

				x.push_front((x[0]*y[1] - x[1]*y[0])/(y[1] - y[0]));
				y.push_front(f(x[0]));

				return x[0];
			}
			X inverse_quadratic(const std::function<X(X)>& f)
			{
				ensure (x.size() >= 3);
				ensure (y[0] != y[1]);
				ensure (y[1] != y[0]);

				X dy01 = y[1] - y[0];
				X dy12 = y[2] - y[1];
				X dy20 = y[0] - y[2];

				x.push_front(x[0]*y[1]*y[2]/(dy01*dy02) + x[1]*y[2]*y[0]/(dy01*dy12) + x[2]*y[0]*y[1]/(dy20*dy12));
				y.push_front(f(x[0]));

				return x[0];
			}
			X newton(const std::function<X(X)>& f, const std::function<X(X)>& df)
			{
				ensure (x.size() >= 1);

				X df0 = df(x[0]);
				ensure (df0 != 0);

				y.push_front(f(x[0]));
				x.push_front(x[0] - y[0]/df);

				return x[0];
			}
		}; 

	} // root1d
} // fms